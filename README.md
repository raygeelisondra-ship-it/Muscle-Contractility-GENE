# Molecular Mechanisms of Muscle Contractility

* **Student Name:** Ray Gee J. Lisondra
* **Course:** Cell and Molecular Biology | Organism: Homo sapiens

---

## Assignment Information
* **Assignment Number:** 5
* **Gene:** MYH3
* **Official Name:** Myosin heavy chain 3
* **Gene Symbol:** MYH3
* **NCBI Gene ID:** 4621
* **Chromosome:** 17p13.1
* **Transcript Accession:** NM_002470.4
* **Variant:** NM_002470.4:c.2015G>A
* **Protein Change:** p.Arg672His
* **ClinVar Variation ID:** 14138
* **Protein Accession:** NP_002461.1
* **Variant Type:** Missense
* **HGVS Change:**
  - Nucleotide: c.2015G>A
  - Protein: p.Arg672His
* **Variant Type & Length:** Single nucleotide variant (substitution)
* **Condition:** Freeman-Sheldon syndrome (Distal arthrogryposis 2A)
* **Aggregate Classification:** Pathogenic
* **Review Status:** 2 gold stars
* **Last Evaluated:** (As recorded in ClinVar database)
* **Date Accessed:** October 27, 2026
* **Submitter:** Clinical laboratories and expert panels via ClinVar
* **Supporting Publications:** PMID: 34664542, 30826400

## Research Question
How does the pathogenic missense variant p.Arg672His in the MYH3 gene structurally and functionally impact embryonic myosin heavy chain motor activity?


## Methods
* **1. Project Setup:**  Install and load required packages, define parameters, and create subfolders
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

* **2. Gene Retrieval:** Search NCBI Gene for human MYH3 and record summary info

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

* **3. Protein Sequence Retrieval:** Fetch transcript record, extract protein accession, and download FASTA

gb <- entrez_fetch(db = "nuccore", id = transcript_accession, rettype = "gb", retmode = "text")
pat <- '/protein_id="NP_[0-9]+[.][0-9]+"' 
fragment <- regmatches(gb, regexpr(pat, gb))
protein_accession <- gsub('^/protein_id="|"$', '', fragment)

protein_fasta <- entrez_fetch(db = "protein", id = protein_accession, rettype = "fasta", retmode = "text")
fasta_file <- file.path("data", paste0(gene_symbol, "_protein.fasta"))
writeLines(protein_fasta, fasta_file)
writeLines(paste("NP accession:", protein_accession), "results/protein_accession.txt")

* **4. Protein Analysis:** Calculate length, count standard amino acids, and generate composition plot

fasta_lines <- readLines(fasta_file)
sequence <- paste(fasta_lines[!grepl("^>", fasta_lines)], collapse = "")
aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))
results <- data.frame(AminoAcid = standard, Count = as.integer(counts))
write.csv(results, "results/amino_acid_composition.csv", row.names = FALSE)

fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") + theme_bw() +
  labs(title = paste(gene_symbol, "amino-acid composition"), x = "Amino acid", y = "Number of residues")
ggsave("figures/01_amino_acid_composition.png", fig1, width = 8, height = 5, dpi = 300)

* **5. Variant Review & Biological Visualization:** Generate linear protein coordinate map for variant position

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

## Results and Figures

* **Transcript Accession:** NM_002470.4
* **Retrieved Protein Accession:** NP_002461.1 (Length: 1,940 amino acids)

### Data Output Files Generated:
* `results/gene_summary.txt`: Contains official gene symbol, NCBI Gene ID (4621), description, and chromosome location (17p13.1).
* `results/protein_accession.txt`: Records the confirmed RefSeq protein accession version.
* `results/amino_acid_composition.csv`: Tabulated counts of all 20 standard amino acids.
* `data/MYH3_protein.fasta`: Reference protein sequence.

### Figures:
1. **Figure 1:** `[Amino Acid Composition Plot](figures/01_amino_acid_composition.png)`
   - **Interpretation:** Displays counts for all 20 standard amino acids in the MYH3 protein, showing that Leucine (L), Glutamic acid (E), and Lysine (K) are the most frequent residues. This profile describes global amino-acid abundance but does not reveal 3D folding or disease mechanisms.

2. **Figure 2:** `[Variant Position Map](figures/02_variant_position.png)`
   - **Interpretation:** Linear coordinate map of the 1,940-amino-acid reference protein highlighting the exact location of residue position 672 where the p.Arg672His missense variant occurs.

3. **Figure 3:** `[Protein Structure and Domain Map](figures/protein_structure_or_splice_diagram.png)`
   - **Interpretation:** UniProt feature viewer and AlphaFold 3D structural model showing the globular motor head domain and structural confidence scores.
  
## Variant interpretation (observed effects, clinical association, and predictions)

The MYH3 c.2015G>A (p.Arg672His) variant is classified as Pathogenic in ClinVar and is strongly associated with Freeman–Sheldon syndrome, also known as distal arthrogryposis type 2A (DA2A). This variant is one of the most commonly reported disease-causing variants in affected individuals and is absent from the general population database, supporting that it is not a benign polymorphism.

The clinical evidence for pathogenicity is strong. Individuals carrying p.Arg672His often show the characteristic features of Freeman–Sheldon syndrome, including a small mouth, pursed lips that create a “whistling-face” appearance, and congenital contractures of the hands and feet. These contractures limit movement from birth and result from abnormal muscle and skeletal development before birth. 
The variant has been repeatedly identified in molecularly confirmed cases of Freeman–Sheldon syndrome, including de novo cases in which it was not present in either parent. This repeated clinical observation links the variant directly to the disorder rather than suggesting a coincidental association.
MYH3 encodes embryonic myosin heavy chain, a motor protein involved in muscle contraction during embryonic and fetal development. Myosin produces force by binding actin filaments and using energy from ATP hydrolysis. The p.Arg672His substitution replaces a conserved positively charged arginine with histidine at position 672 in the globular motor head domain. 
Arg672 lies near the ATP-binding region and interacts with structural elements that move during the myosin cross-bridge cycle. Replacing arginine with histidine may therefore disturb ATP binding, conformational changes, and efficient force production.
Functional studies provide direct evidence that this variant disrupts myosin function. ClinVar reports that functional analyses of the variant protein show defects in myosin function.  Related studies of Freeman–Sheldon-associated myosin mutations, including substitutions at Arg672, show that these residues are close to the ATP-binding pocket and are important for ATPase activity, actin binding, and contractile function. 
Studies of related Arg672 mutations in embryonic myosin have also shown reduced muscle force, prolonged relaxation, and impaired cross-bridge cycling. Although these findings come from related variants and model systems, they support the conclusion that changes at this conserved position can impair muscle contraction. 
Computational tools also predict that p.Arg672His is deleterious. For example, ClinVar reports a REVEL score of 0.814, which indicates a high predicted likelihood of a deleterious effect.  Earlier structure-function analyses likewise predicted that mutations at Arg672 could interfere with ATP binding and myosin catalytic activity. 
These predictions are consistent with the known structure of the myosin motor domain, but they are supportive evidence rather than proof of pathogenicity on their own.
The pathogenic classification of MYH3 p.Arg672His is supported by converging evidence. Clinically, the variant is recurrently found in individuals with Freeman–Sheldon syndrome. Experimentally, it has been shown to cause defects in myosin function. Computationally, it is predicted to be deleterious and may disrupt ATP binding and motor-domain activity. Together, these findings indicate that p.Arg672His impairs embryonic myosin function, reduces effective muscle contraction during development, and contributes to the congenital contractures and craniofacial features of Freeman–Sheldon syndrome. 

### Molecular Causal Chain
Gene (MYH3) ➔ Protein (Embryonic Myosin) ➔ Molecular Function (ATP-Driven Actin Sliding) ➔ Variant (p.Arg672His) ➔ Effect (Impaired Force Generation and Abnormal Contractility)

## Limitations and Validation Experiment

* **Limitations:** Sequence analysis and linear position mapping only provide structural context; they cannot directly measure in vivo motor kinetics, actin-sliding dynamics, or clinical severity variations.
* **Validation Experiment:** An *in vitro* motility assay or recombinant ATPase activity assay comparing wild-type and p.Arg672His mutant myosin proteins to directly measure force generation and enzymatic kinetics.

 * **1. What does your assigned gene encode, and how does this protein contribute to muscle function?**
MYH3 encodes the embryonic heavy chain isoform of myosin (Myosin-3), a molecular motor that powers muscle contraction during fetal development by hydrolyzing ATP and interacting with actin thin filaments.

* **2. What was the exact protein accession and length retrieved with R?**
Accession NP_002461.1 with a total length of 1,940 amino acids.

* **3. What did the amino-acid composition chart show? What can it NOT tell you?**
It showed quantitative distributions of all 20 standard amino acids, with Leucine, Glutamic acid, and Alanine being most frequent[. It cannot determine 3D folding, active site locations, or the functional consequences of missense variants.

* **4. What is your fixed ClinVar variant, and what is its current classification and evidence strength?**
Variant p.Arg672His (c.2015G>A), ClinVar Variation ID 14138, classified as Pathogenic with a 2-star review status, associated with Freeman-Sheldon syndrome.

* **5.What is scientifically supported mechanism linking this variant to muscle function, and what remains uncertain?**
The substitution of a positively charged Arginine with Histidine in the conserved globular motor head domain disrupts electrostatic interactions critical for ATP hydrolysis and cross-bridge cycling. The exact downstream compensatory signaling pathways and variable clinical expressivity remain partially uncertain.

* **6. What additional experiment would help test the proposed mechanism?**
An in vitro actin-sliding motility assay or recombinant ATPase activity assay comparing wild-type and p.Arg672His mutant myosin proteins to measure force generation and enzymatic kinetics.

## (g) Database Links, Accession Versions, Access Dates, and Citations

* **NCBI Gene:** [MYH3 Gene Record (ID: 4621)](https://www.ncbi.nlm.nih.gov/gene/4621) — Accessed October 27, 2026
* **ClinVar Variant:** [ClinVar Variation 14138](https://www.ncbi.nlm.nih.gov/clinvar/variation/14138/) — Accessed October 27, 2026
* **UniProt Entry:** [P11055 (MYH3_HUMAN)](https://www.uniprot.org/uniprotkb/P11055/entry) — Accessed October 27, 2026
* **Supporting Publications:** 
  - Das, S., Kumar, P., Verma, A., Maiti, T. K., & Mathew, S. J. (2019). Myosin heavy chain mutations that cause Freeman-Sheldon syndrome lead to muscle structural and functional defects in Drosophila. *Developmental Biology*, 449(2), 90–98. [PMID: 30826400](https://pubmed.ncbi.nlm.nih.gov/30826400/)
  - Bowman, S., Noble, G., Rahmani, B., Mets, M., Ralay Ranaivo, H., & Castelluccio, V. (2022). A case of blepharophimosis: Freeman Sheldon syndrome. *Ophthalmic Genetics*, 43(1), 130–133. [PMID: 34664542](https://pubmed.ncbi.nlm.nih.gov/34664542/)
* **Reference Citation:** NCBI RefSeq transcript NM_002470.4 and protein NP_002461.1 records.
  
  

  
