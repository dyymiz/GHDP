
# SOURCE FILES in this directory except main.R
sourceDir <- function(path, trace = TRUE, except="main.R") {
  for (nm in setdiff(list.files(path, pattern = "\\.[RrSsQq]$"),except)) {
    if(trace) cat(nm,":")
    source(file.path(path, nm))
    if(trace) cat("\n")
  }
}
#
