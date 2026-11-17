# Run from the repository root. Static inspection only; never sources analyses.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else tempfile('paracocco-validation-')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
files <- list.files(c('code', 'functions', 'data', 'scripts'), pattern = '[.][rR]$', recursive = TRUE, full.names = TRUE)
rows <- list(); calls <- list(); deps <- character()
walk <- function(x, f) {
  if (is.call(x)) {
    nm <- paste(deparse(x[[1]]), collapse = '')
    if (nm %in% c('library', 'require', 'requireNamespace') && length(x) > 1) {
      deps <<- c(deps, as.character(x[[2]]))
    }
    if (is.call(x[[1]]) && as.character(x[[1]][[1]]) %in% c('::', ':::')) deps <<- c(deps, as.character(x[[1]][[2]]))
    if (nm %in% c('source','readRDS','read.csv','readr::read_tsv','readr::read_csv','saveRDS','write.csv','ggsave','ggplot2::ggsave','tiff','pdf','grDevices::tiff','set.seed','setwd','openxlsx::write.xlsx','Matrix::writeMM','writeLines')) {
      calls[[length(calls)+1L]] <<- data.frame(script=f, call=nm, expression=paste(deparse(x, width.cutoff=500), collapse=' '))
    }
    for (i in seq_along(x)) {
      if (!identical(x[[i]], quote(expr = ))) walk(x[[i]], f)
    }
  } else if (is.expression(x) || is.pairlist(x)) for (i in seq_along(x)) {
    if (!identical(x[[i]], quote(expr = ))) walk(x[[i]], f)
  }
}
for (f in files) {
  err <- NULL
  ex <- tryCatch(parse(f), error=function(e) {err <<- conditionMessage(e); NULL})
  rows[[length(rows)+1L]] <- data.frame(script=f, syntax=if(is.null(err)) 'STATIC PASS' else 'REVIEW REQUIRED', detail=if(is.null(err)) 'Parsed; runtime not implied' else err)
  if (!is.null(ex)) walk(ex, f)
}
write.table(do.call(rbind, rows), file.path(out,'r_syntax.tsv'), sep='\t', quote=TRUE, row.names=FALSE)
write.table(do.call(rbind, calls), file.path(out,'r_calls.tsv'), sep='\t', quote=TRUE, row.names=FALSE)
ip <- installed.packages()
deps <- sort(unique(c(deps, 'DDRTree','igraph','mgcv','Biobase')))
env <- data.frame(package=deps, installed=deps %in% rownames(ip), version=vapply(deps, function(p) if(p %in% rownames(ip)) ip[p,'Version'] else NA_character_, character(1)))
write.table(env, file.path(out,'r_environment.tsv'), sep='\t', quote=TRUE, row.names=FALSE)
capture.output(sessionInfo(), file=file.path(out,'sessionInfo.txt'))
cat('Static audit written to',out,'\n')
if (any(vapply(rows, function(x) x$syntax != 'STATIC PASS', logical(1)))) quit(status=1)
