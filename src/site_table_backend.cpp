#include <Rcpp.h>
#include "site_table.h"

using namespace Rcpp;
using riskmap_site::Evaluator;
using riskmap_site::Moments;
using riskmap_site::SiteTable;
using riskmap_site::Unit;

// Tables of records with their independent effects integrated out
// (site_table.h), for the adaptive MCML engine with integrate_noise. A table
// object is created once per engine, for the distinct pairs (y, m) of the
// records and of the units of its nested patterns, and holds the record tables
// of up to 32 variances and the nested tables of up to 16 pairs of variances,
// and for groups of units whose linear predictors differ, such as the
// children of a household, the tables of groups and of locations of groups at
// the current state of their relative linear predictors.

namespace {

SiteTable& site_table_from(SEXP table) {
  XPtr<SiteTable> pointer(table);
  if (pointer.get() == nullptr) stop("The table of integrated sites is empty.");
  return *pointer;
}

// The engine's records from a matrix with one row per record and the columns
// kind (0, 1 or 2), index, a and b of riskmap_site::Unit, the last three
// 1-based.
std::vector<Unit> units_from(const IntegerMatrix& units) {
  if (units.ncol() != 4) stop("The records need four columns.");
  std::vector<Unit> out(units.nrow());
  for (int r = 0; r < units.nrow(); ++r) {
    out[r] = {units(r, 0), units(r, 1) - 1, units(r, 2) - 1, units(r, 3) - 1};
  }
  return out;
}

}  // namespace

// `z` and `log_w` are the nodes and log weights of Gauss-Hermite quadrature
// for exp(-z^2). `custom` is NULL for the canonical links, or for a custom
// link a list with the lattice `low` and `step` and the arrays A, A1, A2, B,
// B1 and B2 of its kernel (site_table.h). `patterns` holds the 1-based pairs
// of the units of each nested pattern, and `summed` the 1-based pair of each
// pattern's units summed into one record. `groups` holds the 1-based pairs of
// each group's members, `locations` the 1-based groups of each location of
// groups, and `whole` the 1-based group of all of each location's members.
// `concave` says that a custom link is log-concave (SiteTable).
// [[Rcpp::export]]
SEXP cpp_site_table_create(const NumericVector& y, const NumericVector& m,
                           const bool poisson, const NumericVector& z,
                           const NumericVector& log_w, SEXP custom,
                           const List& patterns, const IntegerVector& summed,
                           const List& groups, const List& locations,
                           const IntegerVector& whole, const bool concave) {
  riskmap_site::Rule rule;
  rule.z.assign(z.begin(), z.end());
  rule.log_w.assign(log_w.begin(), log_w.end());
  riskmap_site::CustomArrays arrays;
  const bool has_custom = !Rf_isNull(custom);
  if (has_custom) {
    const List parts(custom);
    arrays.low = as<double>(parts["low"]);
    arrays.step = as<double>(parts["step"]);
    const char* names[] = {"A", "A1", "A2", "B", "B1", "B2"};
    std::vector<double>* targets[] = {&arrays.A, &arrays.A1, &arrays.A2,
                                      &arrays.B, &arrays.B1, &arrays.B2};
    for (int i = 0; i < 6; ++i) {
      const NumericVector values = parts[names[i]];
      targets[i]->assign(values.begin(), values.end());
      if (targets[i]->size() != arrays.A.size() || values.size() < 2) {
        stop("The arrays of a custom link must have one common length.");
      }
    }
  }
  if (patterns.size() != summed.size()) {
    stop("Each pattern needs the pair of its summed units.");
  }
  std::vector<std::vector<int>> pattern_pairs(patterns.size());
  std::vector<int> summed_pairs(summed.size());
  for (R_xlen_t p = 0; p < patterns.size(); ++p) {
    const IntegerVector pairs = patterns[p];
    for (const int pair : pairs) {
      if (pair < 1 || pair > y.size()) stop("A pattern names a missing pair.");
      pattern_pairs[p].push_back(pair - 1);
    }
    if (summed[p] < 1 || summed[p] > y.size()) {
      stop("A pattern names a missing summed pair.");
    }
    summed_pairs[p] = summed[p] - 1;
  }
  std::vector<std::vector<int>> group_pairs(groups.size());
  for (R_xlen_t g = 0; g < groups.size(); ++g) {
    const IntegerVector pairs = groups[g];
    if (pairs.size() == 0) stop("A group has no members.");
    for (const int pair : pairs) {
      if (pair < 1 || pair > y.size()) stop("A group names a missing pair.");
      group_pairs[g].push_back(pair - 1);
    }
  }
  if (locations.size() != whole.size()) {
    stop("Each location of groups needs the group of all its members.");
  }
  std::vector<std::vector<int>> location_groups(locations.size());
  std::vector<int> whole_groups(whole.size());
  for (R_xlen_t l = 0; l < locations.size(); ++l) {
    const IntegerVector members = locations[l];
    if (members.size() == 0) stop("A location of groups has no groups.");
    for (const int g : members) {
      if (g < 1 || g > groups.size()) stop("A location names a missing group.");
      location_groups[l].push_back(g - 1);
    }
    if (whole[l] < 1 || whole[l] > groups.size()) {
      stop("A location names a missing group of all its members.");
    }
    whole_groups[l] = whole[l] - 1;
  }
  SiteTable* table = new SiteTable(std::vector<double>(y.begin(), y.end()),
                                   std::vector<double>(m.begin(), m.end()),
                                   poisson, rule, std::move(arrays),
                                   has_custom, std::move(pattern_pairs),
                                   std::move(summed_pairs), group_pairs,
                                   std::move(location_groups),
                                   std::move(whole_groups), concave);
  return XPtr<SiteTable>(table, true);
}

// Sets the linear predictors of the groups' members, in the order of their
// groups, which the group and location integrals resolved from then on take:
// first each member's relative to its record, then each member's relative to
// the first member of its group (SiteTable::set_relative()).
// [[Rcpp::export]]
void cpp_site_set_relative(SEXP table, const NumericVector& relative) {
  SiteTable& sites = site_table_from(table);
  sites.set_relative(std::vector<double>(relative.begin(), relative.end()));
}

// The log-likelihood of each record and, for `order` 1 or 2, its derivatives,
// for a matrix of linear predictors, at the variances `classes` of the
// records' classes. `units` describes one record for each row of `eta`
// (`by_row`) or each column.
// [[Rcpp::export]]
List cpp_site_evaluate(SEXP table, const NumericVector& classes,
                       const IntegerMatrix& units, const NumericMatrix& eta,
                       const bool by_row, const int order) {
  SiteTable& sites = site_table_from(table);
  const int rows = eta.nrow(), columns = eta.ncol();
  if (units.nrow() != (by_row ? rows : columns)) {
    stop("One record description is needed for each record.");
  }
  const std::vector<Unit> described = units_from(units);
  Evaluator evaluator(sites, std::vector<double>(classes.begin(), classes.end()),
                      described);
  NumericMatrix h(rows, columns), h1(rows, columns), h2(rows, columns);
  for (int column = 0; column < columns; ++column) {
    for (int row = 0; row < rows; ++row) {
      const Moments out = evaluator.moments(by_row ? row : column,
                                            eta(row, column), order);
      h(row, column) = out.h;
      h1(row, column) = out.h1;
      h2(row, column) = out.h2;
    }
  }
  if (order == 0) return List::create(Named("h") = h);
  return List::create(Named("h") = h, Named("h1") = h1, Named("h2") = h2);
}

// The sum over the columns (records) of the log-likelihood, for each row
// (draw) of `eta`.
// [[Rcpp::export]]
NumericVector cpp_site_row_sums(SEXP table, const NumericVector& classes,
                                const IntegerMatrix& units,
                                const NumericMatrix& eta) {
  SiteTable& sites = site_table_from(table);
  const int rows = eta.nrow(), columns = eta.ncol();
  if (units.nrow() != columns) stop("One record description is needed for each record.");
  const std::vector<Unit> described = units_from(units);
  Evaluator evaluator(sites, std::vector<double>(classes.begin(), classes.end()),
                      described);
  NumericVector out(rows);
  for (int column = 0; column < columns; ++column) {
    for (int row = 0; row < rows; ++row) {
      out[row] += evaluator.log_likelihood(column, eta(row, column));
    }
  }
  return out;
}

// For each draw (row of `eta`) and record (column), the derivative of a
// nested pattern's or a location of groups' log-likelihood in its inner
// variance, and zero for the other records.
// [[Rcpp::export]]
NumericMatrix cpp_site_inner_slopes(SEXP table, const NumericVector& classes,
                                    const IntegerMatrix& units,
                                    const NumericMatrix& eta) {
  SiteTable& sites = site_table_from(table);
  const int rows = eta.nrow(), columns = eta.ncol();
  if (units.nrow() != columns) stop("One record description is needed for each record.");
  const std::vector<Unit> described = units_from(units);
  Evaluator evaluator(sites, std::vector<double>(classes.begin(), classes.end()),
                      described);
  NumericMatrix out(rows, columns);
  for (int column = 0; column < columns; ++column) {
    if (described[column].kind != 2 && described[column].kind != 4) continue;
    for (int row = 0; row < rows; ++row) {
      out(row, column) = evaluator.inner_slope(column, eta(row, column));
    }
  }
  return out;
}

// The number of pairs, patterns, groups and locations of groups, and of
// tables built and cells filled so far, for records, nested patterns, groups
// and locations of groups.
// [[Rcpp::export]]
List cpp_site_table_info(SEXP table) {
  SiteTable& sites = site_table_from(table);
  return List::create(
    Named("pairs") = sites.pairs(),
    Named("tables_built") = sites.tables_built(),
    Named("cells_filled") = static_cast<double>(sites.cells_filled()),
    Named("patterns") = sites.pattern_count(),
    Named("nested_tables_built") = sites.nested_tables_built(),
    Named("nested_cells_filled") = static_cast<double>(sites.nested_cells_filled()),
    Named("groups") = sites.group_count(),
    Named("group_locations") = sites.location_count(),
    Named("group_tables_built") = sites.group_tables_built(),
    Named("group_cells_filled") = static_cast<double>(sites.group_cells_filled()),
    Named("location_tables_built") = sites.location_tables_built(),
    Named("location_cells_filled") = static_cast<double>(sites.location_cells_filled()));
}
