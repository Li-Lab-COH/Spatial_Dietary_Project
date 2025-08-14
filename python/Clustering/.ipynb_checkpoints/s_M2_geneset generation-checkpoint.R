# ------------------------------------------------------------
# 0) Install / load needed packages
# ------------------------------------------------------------
if (!requireNamespace("Seurat",   quietly=TRUE)) install.packages("Seurat")
if (!requireNamespace("biomaRt",  quietly=TRUE)) install.packages("biomaRt")
if (!requireNamespace("jsonlite", quietly=TRUE)) install.packages("jsonlite")

library(Seurat)
library(biomaRt)
library(jsonlite)

# ------------------------------------------------------------
# 1) Pull in the 2019-updated human CC gene sets from Seurat
# ------------------------------------------------------------
# (this object ships with Seurat v4+)
s.human   <- cc.genes.updated.2019$s.genes
g2m.human <- cc.genes.updated.2019$g2m.genes

# ------------------------------------------------------------
# 2) A robust converter that tries multiple Ensembl mirrors
# ------------------------------------------------------------
ensembl_hosts <- c(
  "https://www.ensembl.org",
  "https://uswest.ensembl.org",
  "https://useast.ensembl.org",
  "https://asia.ensembl.org"
)

convertHumanToMouse <- function(human_syms, hosts=ensembl_hosts) {
  for (h in hosts) {
    message("Trying Ensembl host: ", h)
    # attempt to connect
    hm <- try(useMart("ensembl", host=h, dataset="hsapiens_gene_ensembl"), silent=TRUE)
    mm <- try(useMart("ensembl", host=h, dataset="mmusculus_gene_ensembl"), silent=TRUE)
    if (inherits(hm, "Mart") && inherits(mm, "Mart")) {
      out <- try(getLDS(
        attributes   = "hgnc_symbol",
        filters      = "hgnc_symbol",
        values       = unique(human_syms),
        mart         = hm,
        attributesL  = "mgi_symbol",
        martL        = mm,
        uniqueRows   = TRUE
      ), silent=TRUE)
      if (!inherits(out, "try-error")) {
        mouse_syms <- unique(out[,2])
        mouse_syms <- mouse_syms[mouse_syms!="" & !is.na(mouse_syms)]
        message("Mapped ", length(mouse_syms), " symbols via ", h)
        return(mouse_syms)
      }
    }
  }
  stop("All Ensembl hosts failed to map human → mouse symbols.")
}

# ------------------------------------------------------------
# 3) Perform the conversion
# ------------------------------------------------------------
s.mouse   <- convertHumanToMouse(s.human)
g2m.mouse <- convertHumanToMouse(g2m.human)

# ------------------------------------------------------------
# 4) Save out to disk for easy Python use
# ------------------------------------------------------------
# 4a) JSON dict
out_list <- list(s_genes   = s.mouse,
                 g2m_genes = g2m.mouse)

write_json(out_list,
           path       = "cc_mouse_cycle_genes.json",
           auto_unbox = TRUE,
           pretty     = TRUE)

# 4b) Plain‑text, one gene per line
writeLines(s.mouse,   "s_genes_mouse.txt")
writeLines(g2m.mouse, "g2m_genes_mouse.txt")

message("✅ Files written:\n",
        " • cc_mouse_cycle_genes.json\n",
        " • s_genes_mouse.txt\n",
        " • g2m_genes_mouse.txt")
