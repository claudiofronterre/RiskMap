#include <Rcpp.h>
#include "site_table.h"

using namespace Rcpp;
using riskmap_site::Moments;
using riskmap_site::SiteTable;
using riskmap_site::Table;

// Tables of records with the nugget integrated out (site_table.h), for the
// adaptive MCML engine with integrate_noise = TRUE. A table object is created
// once per engine, for the distinct pairs (y, m) of the records whose nugget
// is integrated, and holds the tables of up to 16 values of the nugget
// variance v.

namespace {

SiteTable& site_table_from(SEXP table) {
  XPtr<SiteTable> pointer(table);
  if (pointer.get() == nullptr) stop("The table of integrated sites is empty.");
  return *pointer;
}

}  // namespace

// `z` and `log_w` are the nodes and log weights of Gauss-Hermite quadrature
// for exp(-z^2). `custom` is NULL for the canonical links, or for a custom
// link a list with the lattice `low` and `step` and the arrays A, A1, A2, B,
// B1 and B2 of its kernel (site_table.h).
// [[Rcpp::export]]
SEXP cpp_site_table_create(const NumericVector& y, const NumericVector& m,
                           const bool poisson, const NumericVector& z,
                           const NumericVector& log_w, SEXP custom) {
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
  SiteTable* table = new SiteTable(std::vector<double>(y.begin(), y.end()),
                                   std::vector<double>(m.begin(), m.end()),
                                   poisson, rule, std::move(arrays),
                                   has_custom);
  return XPtr<SiteTable>(table, true);
}

// h(eta) and, for `order` 1 or 2, its derivatives, for a matrix of linear
// predictors. `pair` holds a 1-based pair for each row (`by_row`) or each
// column, negative for the plain kernel of that pair. With v = 0 every record
// takes the plain kernel, as for a nugget held at zero.
// [[Rcpp::export]]
List cpp_site_evaluate(SEXP table, const double v, const IntegerVector& pair,
                       const NumericMatrix& eta, const bool by_row,
                       const int order) {
  SiteTable& sites = site_table_from(table);
  const int rows = eta.nrow(), columns = eta.ncol();
  if (pair.size() != (by_row ? rows : columns)) {
    stop("One pair is needed for each record.");
  }
  Table* current = v > 0 ? &sites.table(v) : nullptr;
  NumericMatrix h(rows, columns), h1(rows, columns), h2(rows, columns);
  for (int column = 0; column < columns; ++column) {
    for (int row = 0; row < rows; ++row) {
      const int index = pair[by_row ? row : column];
      const Moments out = (current == nullptr || index < 0) ?
        sites.plain(std::abs(index) - 1, eta(row, column)) :
        sites.evaluate(*current, index - 1, eta(row, column), order);
      h(row, column) = out.h;
      h1(row, column) = out.h1;
      h2(row, column) = out.h2;
    }
  }
  if (order == 0) return List::create(Named("h") = h);
  return List::create(Named("h") = h, Named("h1") = h1, Named("h2") = h2);
}

// The sum over the columns (records) of h, for each row (draw) of `eta`.
// [[Rcpp::export]]
NumericVector cpp_site_row_sums(SEXP table, const double v,
                                const IntegerVector& pair,
                                const NumericMatrix& eta) {
  SiteTable& sites = site_table_from(table);
  const int rows = eta.nrow(), columns = eta.ncol();
  if (pair.size() != columns) stop("One pair is needed for each record.");
  Table* current = v > 0 ? &sites.table(v) : nullptr;
  NumericVector out(rows);
  for (int column = 0; column < columns; ++column) {
    for (int row = 0; row < rows; ++row) {
      out[row] += riskmap_site::record_log_likelihood(
        sites, current, pair[column], eta(row, column)
      );
    }
  }
  return out;
}

// The number of pairs, of tables built and of cells filled so far.
// [[Rcpp::export]]
List cpp_site_table_info(SEXP table) {
  SiteTable& sites = site_table_from(table);
  return List::create(Named("pairs") = sites.pairs(),
                      Named("tables_built") = sites.tables_built(),
                      Named("cells_filled") = static_cast<double>(sites.cells_filled()));
}
