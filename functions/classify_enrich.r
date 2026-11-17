# Shared helper: GO program signatures and annotation helpers; annotate_th17_terms is explicit.
# Inputs: caller arguments/session objects; no automatic data loading.
# Outputs: functions and documented constants; plots only when called.
# Dependencies/prerequisites: see docs/SCRIPT_INVENTORY.md and SOFTWARE_REQUIREMENTS.md.




get_ordering_genes <- function(cds) {

  fd <- fData(cds)

  stopifnot("use_for_ordering" %in% colnames(fd))

  unique(
    rownames(fd)[which(fd$use_for_ordering)]
  )
}




count_signature_genes <- function(
    gene_string,
    signature
) {

  if (
    is.na(gene_string) ||
    !nzchar(gene_string)
  ) {
    return(0L)
  }

  term_genes <- unique(
    strsplit(
      gene_string,
      split = "/",
      fixed = TRUE
    )[[1]]
  )

  length(
    intersect(
      term_genes,
      signature
    )
  )
}

get_signature_genes <- function(
    gene_string,
    signature
) {

  if (
    is.na(gene_string) ||
    !nzchar(gene_string)
  ) {
    return("")
  }

  term_genes <- unique(
    strsplit(
      gene_string,
      split = "/",
      fixed = TRUE
    )[[1]]
  )

  paste(
    intersect(
      term_genes,
      signature
    ),
    collapse = "/"
  )
}

# Optional annotation, called only after go_table and th17_signature exist.
annotate_th17_terms <- function(go_table, th17_signature) {
go_table <- go_table |>
  dplyr::mutate(

    N_Th17_genes = vapply(
      .data$geneID,
      count_signature_genes,
      integer(1),
      signature = th17_signature
    ),

    Th17_genes = vapply(
      .data$geneID,
      get_signature_genes,
      character(1),
      signature = th17_signature
    ),

    has_Il17a = stringr::str_detect(
      paste0("/", .data$geneID, "/"),
      stringr::fixed("/Il17a/")
    )
  )
  go_table
}






signatures <- list(

  Less_differentiated = c(
    "Lef1",
    "Klf3",
    "Sell",
    "Foxp1",
    "S1pr1",
    "Klf2",
    "Ccr7",
    "Tob1",
    "Tcf7",
    "Tob2",
    "Zfp36l2",
    "Btg2"
  ),

  T_cell_activation = c(
    "Tnfrsf1b",
    "Icos",
    "Cyld",
    "Nr4a1",
    "Itk",
    "Tnfaip3",
    "Tnfsf8",
    "Cd28",
    "Zap70",
    "Lat"
  ),

  Migration_and_adhesion = c(
    "Rgs1",
    "Cxcr4",
    "Ccr2",
    "Ccr7",
    "S1pr1",
    "Spn",
    "Adam8",
    "Itgb1",
    "Itga4",
    "S1pr4",
    "Rac2",
    "Ripor2",
    "Cxcr6",
    "S100a4",
    "Ccl5"
  ),

  Th17_associated = c(
    "Il17a",
    "Maf",
    "Hif1a",
    "Bhlhe40"
  ),

  Treg_associated = c(
    "Ctla4",
    "Il2ra",
    "Tmem176a",
    "Tmem176b",
    "Tgfbr2"
  )
)
















split_gene_ids <- function(gene_string) {

  if (
    length(gene_string) == 0 ||
    is.na(gene_string) ||
    !nzchar(gene_string)
  ) {
    return(character(0))
  }

  genes <- unlist(
    strsplit(
      gene_string,
      split = "/",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  genes <- trimws(genes)

  unique(
    genes[
      !is.na(genes) &
        nzchar(genes)
    ]
  )
}


count_signature_genes <- function(
    gene_string,
    signature
) {

  term_genes <- split_gene_ids(gene_string)

  length(
    intersect(
      term_genes,
      signature
    )
  )
}


get_signature_genes <- function(
    gene_string,
    signature
) {

  term_genes <- split_gene_ids(gene_string)

  matched_genes <- intersect(
    term_genes,
    signature
  )

  paste(
    matched_genes,
    collapse = "/"
  )
}


contains_gene <- function(
    gene_string,
    gene
) {

  gene %in% split_gene_ids(gene_string)
}
