.onLoad <- function(libname, pkgname) {
  if (
    requireNamespace("prolfqua", quietly = TRUE) &&
      exists("register_facade", where = asNamespace("prolfqua"), inherits = FALSE)
  ) {
    prolfqua::register_facade(
      "saint",
      class = "ContrastsSAINTFacade",
      needs = "aggregated",
      package = "prolfquasaint",
      needs_saint_annotation = TRUE
    )
  }
  invisible()
}
