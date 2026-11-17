# Execute a reviewed output-only panel in a fresh temporary directory.
# Rscript --vanilla scripts/smoke_test_panel.R code/Figure_3/figure_3I/code.r
args <- commandArgs(trailingOnly = TRUE)
allowed <- c('code/Figure_1/figure_1B/code.r', 'code/Figure_1/figure_1C/code.r',
 'code/Figure_1/figure_1D/code.r', 'code/Figure_1/Figure_1E_and_1F/code.r',
 'code/Figure_1/Figure_1G/code.r', 'code/Figure_2/figure_2b/code.r',
 'code/Figure_2/figure_2D-E/code.r', 'code/Figure_2/figure_2F/code.r',
 'code/Figure_3/figure_3I/code.r', 'code/Figure_3/figure_3J/code.r',
 'code/Figure_3/figure_3H/heatmap.r', 'code/Figure_3/figure_3D-G/relative_effect.r')
if (length(args) != 1L || !args[[1]] %in% allowed) stop('Choose one reviewed output-only panel from the allowlist in this script.')
root <- normalizePath('.')
panel <- args[[1]]
out <- tempfile('paracocco-panel-'); dir.create(out)
stopifnot(file.symlink(file.path(root,'data'),file.path(out,'data')),
          file.symlink(file.path(root,'functions'),file.path(out,'functions')))
setwd(out)
cat('Output directory:',out,'\nSTART',panel,'\n')
tryCatch({
  source(file.path(root,panel), local = .GlobalEnv, print.eval = TRUE)
  if (panel == 'code/Figure_3/figure_3H/heatmap.r') {
    source(file.path(root,'code/Figure_3/figure_3H/alluvial.r'), local = .GlobalEnv, print.eval = TRUE)
    cat('PASS code/Figure_3/figure_3H/alluvial.r\n')
  }
  cat('PASS',panel,'\n')
}, error=function(e) {cat('BLOCKED',panel,conditionMessage(e),'\n'); quit(status=1)})
graphics.off()
print(file.info(list.files(out,recursive=TRUE,include.dirs=FALSE))[,'size',drop=FALSE])
