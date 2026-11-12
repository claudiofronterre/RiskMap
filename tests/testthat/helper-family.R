# The family object of a family string with a link, or with the kernel of
# integrated records (make_integrated_kernel()), the two forms that the engine's
# functions took before step 1d of #253. An object is returned as it is.
family_object <- function(family, linkf = NULL) {
  if (inherits(family, "riskmap_family")) return(family)
  if (is_integrated_link(linkf)) return(make_family(family, kernel = linkf))
  make_family(family, linkf)
}
