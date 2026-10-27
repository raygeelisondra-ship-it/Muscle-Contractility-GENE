# Molecular Mechanisms of Muscle Contractility

* **Name:** Ray Gee J. Lisondra
* **Assignment Number:** 05
* **Assigned Gene:** MYH3 (Myosin Heavy Chain 3)
* **Transcript:** NM_002470.4
* **ClinVar ID:** 14138
* **Associated Disease:** Freeman-Sheldon syndrome (Distal arthrogryposis type 2A)

## Assigned Gene and Disease
* **Gene Symbol:** MYH3
* **NCBI Gene ID:** 4621
* **Full Gene Name:** Myosin heavy chain 3
* **Chromosome:** Chromosome 17 (17p13.1)
* **Normal Role in Muscle:** Acts as a molecular motor using ATP hydrolysis to drive actin-myosin sliding and force generation during fetal muscle development.

## Research Question
How does the pathogenic missense variant p.Arg672His in the MYH3 gene structurally and functionally impact embryonic myosin heavy chain motor activity?

## 3. R Analysis Workflow (`analysis.R`)
```r
# Run installation only once on your computer:
install.packages(c("rentrez", "ggplot2"))
library(rentrez)
library(ggplot2)

gene_symbol <- "MYH3"
transcript_accession <- "NM_002470.4" 
clinvar_id <- "14138" 
variant_protein <- "p.Arg672His"

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}

# Retrieve NCBI Gene Record
hits <- entrez_search(
  db = "gene", term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"), retmax = 20
)
if (length(hits$ids) == 0) stop("Gene not found")

gene_id <- NA_character_
for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id; gene_info <- rec; break
  }
}
if (is.na(gene_id)) stop("Exact human gene not found")

writeLines(c(
  paste("Gene symbol:", gene_info$name), 
  paste("NCBI Gene ID:", gene_id), 
  paste("Official name:", gene_info$description), 
  paste("Chromosome:", gene_info$chromosome)
), "results/gene_summary.txt")

# Retrieve Protein Sequence from Transcript
gb <- entrez_fetch(db = "nuccore", id = transcript_accession, rettype = "gb", retmode = "text")
pat <- '/protein_id="NP_[0-9]+[.][0-9]+"' 
fragment <- regmatches(gb, regexpr(pat, gb))
protein_accession <- gsub('^/protein_id="|"$', '', fragment)

protein_fasta <- entrez_fetch(db = "protein", id = protein_accession, rettype = "fasta", retmode = "text")
fasta_file <- file.path("data", paste0(gene_symbol, "_protein.fasta"))
writeLines(protein_fasta, fasta_file)
writeLines(paste("NP accession:", protein_accession), "results/protein_accession.txt")

# Amino-Acid Composition Analysis
fasta_lines <- readLines(fasta_file)
sequence <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))
results <- data.frame(AminoAcid = standard, Count = as.integer(counts))
write.csv(results, "results/amino_acid_composition.csv", row.names = FALSE)

# Generate Amino-Acid Composition Bar Plot
fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") + theme_bw() +
  labs(title = paste(gene_symbol, "amino-acid composition"), x = "Amino acid", y = "Number of residues")
ggsave("figures/01_amino_acid_composition.png", fig1, width = 8, height = 5, dpi = 300)

# Generate Variant Position Map
if (nzchar(variant_protein)) {
  matched <- regmatches(variant_protein, regexpr("[0-9]+", variant_protein))
  position <- as.integer(matched)
  if (!is.na(position) && position >= 1 && position <= protein_length) {
    fig2 <- ggplot() +
      geom_segment(aes(x = 1, xend = protein_length, y = 0, yend = 0), linewidth = 2) +
      geom_point(aes(x = position, y = 0), size = 4) +
      annotate("text", x = position, y = 0.25, label = variant_protein, size = 4) +
      xlim(0, protein_length) + ylim(-0.6, 0.6) +
      theme_void() +
      labs(title = paste(gene_symbol, "reference protein"), x = "Amino-acid position") +
      theme(axis.title.x = element_text())
    ggsave("figures/02_variant_position.png", fig2, width = 8, height = 2.5, dpi = 300)
  }
}


  

  
