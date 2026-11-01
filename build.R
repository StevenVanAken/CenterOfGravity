required <- c(
  "shinylive",
  "httpuv",
  "shiny",
  "leaflet",
  "DT",
  "geosphere",
  "lpSolve"
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
if (!file.exists("app/app.R")) {
  stop(
    "Cannot find app/app.R. Current working directory: ",
    getwd(),
    call. = FALSE
  )
}

shinylive::export(
  appdir = "app",
  destdir = "docs",
  wasm_packages = TRUE
)

# Tell GitHub Pages to serve the export without Jekyll processing.
file.create("docs/.nojekyll")

# Confirm the website entry point exists.
stopifnot(file.exists("docs/index.html"))

message("Export ready. Preview with:")
message('httpuv::runStaticServer("docs")')
