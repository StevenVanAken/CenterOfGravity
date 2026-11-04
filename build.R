# Run from the warehouse-cog project root.
setwd("C:/Users/Steve/OneDrive/R/warehouse-cog/")

required <- c(
  "shinylive",
  "httpuv",
  "shiny",
  "leaflet",
  "DT",
  "geosphere",
  "lpSolve",
  "xfun"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing)) {
  stop(
    "Install these packages before building: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

if (!file.exists("C:/Users/Steve/OneDrive/R/warehouse-cog/app/app.R")) {
  stop(
    "Expected app/app.R. Run this script from the project root.",
    call. = FALSE
  )
}

# Export application, web runtime, and available WASM package binaries.
shinylive::export(
  appdir = "app",
  destdir = "site",
  wasm_packages = TRUE
)

message("Export finished. Preview with:")
message('httpuv::runStaticServer("site")')