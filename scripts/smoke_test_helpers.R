# Source shared helpers and apply compatibility patches twice; does not run analyses.
suppressPackageStartupMessages({library(Seurat);library(monocle);library(dplyr);library(ggplot2)})
for(f in list.files('functions',pattern='[.]r$',recursive=TRUE,full.names=TRUE)) {
 cat('START',f,'\n')
 tryCatch({source(f,local=new.env(parent=globalenv())); cat('PASS helper source',f,'\n')},error=function(e)cat('BLOCKED',f,conditionMessage(e),'\n'))
}
source('functions/patch/patch.r')
tryCatch({patch_monocle2(); patch_monocle2();cat('PASS patch application and repeated application\n')},error=function(e)cat('BLOCKED patch',conditionMessage(e),'\n'))
