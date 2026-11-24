library(data.table)

aln_block_length=1000
aln_qual=40

aln <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/candidate_regions/candidate_seltot_sum.paf", fill=TRUE)
head(aln)

aln<-aln[aln$V11 > aln_block_length & aln$V12 >= aln_qual,]
head(aln)

regions<-paste0(aln$V6,":",as.integer(aln$V8),":",as.integer(aln$V9))

gff_lines <- readLines("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/candidate_regions/GCF_000633615.1_Guppy_female_1.0_MT_genomic-2.gff")
gff_lines <- gff_lines[!startsWith(gff_lines, "#")]

gff <- fread(text = gff_lines,  sep = "\t",  header = FALSE,  quote = "",
  col.names = c(  "seqid", "source", "type", "start", "end",  "score", "strand", "phase", "attributes"))


# Convert to data.table
aln <- as.data.table(aln)
gff <- as.data.table(gff)

# Give each alignment a unique ID
aln[, alignment_id := .I]

# PAF reference coordinates are 0-based, half-open
aln[, `:=`(
  region_start = as.integer(V8) + 1L,
  region_end   = as.integer(V9),
  reference    = V6
)]

# Extract gene features from GFF
genes <- gff[
  type == "gene",
  .(
    reference = seqid,
    gene_start = start,
    gene_end = end,
    attributes
  )
]

# Find overlaps using the original query sequence name
regions_dt <- aln[, .(
  query_sequence = V1,
  reference,
  region_start,
  region_end
)]

setkey(regions_dt, reference, region_start, region_end)
setkey(genes, reference, gene_start, gene_end)

gene_hits <- foverlaps(
  regions_dt,
  genes,
  by.x = c("reference", "region_start", "region_end"),
  by.y = c("reference", "gene_start", "gene_end"),
  type = "any",
  nomatch = 0L
)

# Keep unique attributes per query sequence
gene_hits <- unique(
  gene_hits[, .(query_sequence, attributes)]
)
head(gene_hits)
nrow(gene_hits)

# Save results
fwrite(gene_hits, "Documents/FIBR_temporal/am_nat_revisions/af_sum_change/candidate_regions/candidate_seltot_sum.annot",  sep = "\t")


#GO analysis
extract_ncbi_id <- function(x) {
  out <- sub(".*GeneID:([0-9]+).*", "\\1", x)
  out[!grepl("GeneID:[0-9]+", x)] <- NA_character_
  out
}

# All genes in the genome
genome_genes <- unique(na.omit(
  extract_ncbi_id(gff[type == "gene"]$attributes)
))
head(genome_genes)

# Genes overlapping candidate regions
candidate_genes <- unique(na.omit(
  extract_ncbi_id(gene_hits$attributes)
))
head(candidate_genes)

# Restrict to genes with GO annotations
gene2go <- read.table("Documents/FIBR_temporal/am_nat_revisions/af_sum_change/candidate_regions/guppy_gene2go.tsv", header=T)
head(gene2go)
genome_annotated <- intersect(
  genome_genes, gene2go$gene_id
)

candidate_annotated <- intersect(
  candidate_genes, genome_annotated
)

cat("Genes in GFF:", length(genome_genes), "\n")
cat("Candidate genes:", length(candidate_genes), "\n")
cat("Genome genes with GO annotations:",
    length(genome_annotated), "\n")
cat("Candidate genes with GO annotations:",
    length(candidate_annotated), "\n")
head(gene_hits)

# -----------------------------------------
# 1. Prepare candidate and background genes
# -----------------------------------------

# Use character IDs consistently
gene2go$gene_id <- as.character(gene2go$gene_id)

candidate_genes <- unique(as.character(candidate_genes))
genome_genes <- unique(as.character(genome_genes))

# Restrict to genes with GO annotations
annotated_genes <- unique(gene2go$gene_id)

background <- intersect(genome_genes, annotated_genes)
candidates <- intersect(candidate_genes, background)

cat("Annotated background genes:", length(background), "\n")
cat("Annotated candidate genes:", length(candidates), "\n")

if (length(candidates) == 0) {
  stop("No candidate genes match the annotated background.")
}


# -----------------------------------------
# 2. Perform GO enrichment
# -----------------------------------------

# Restrict annotations to the genomic background
go_map <- unique(
  gene2go[gene2go$gene_id %in% background,
          c("gene_id", "go_id")]
)

# Genes annotated to each GO term
go_sets <- split(go_map$gene_id, go_map$go_id)
go_sets <- lapply(go_sets, unique)

N <- length(background)
n <- length(candidates)

results <- lapply(names(go_sets), function(go) {
  
  genes_with_term <- go_sets[[go]]
  
  # Number of background genes annotated to this term
  K <- length(intersect(genes_with_term, background))
  
  # Number of candidate genes annotated to this term
  overlap <- intersect(genes_with_term, candidates)
  k <- length(overlap)
  
  # Enrichment p-value
  p <- phyper(
    q = k - 1,
    m = K,
    n = N - K,
    k = n,
    lower.tail = FALSE
  )
  
  data.frame(
    GO_ID = go,
    background_genes = K,
    candidate_genes = k,
    p_value = p,
    fold_enrichment = (k / n) / (K / N),
    overlapping_gene_ids = paste(overlap, collapse = ";"),
    stringsAsFactors = FALSE
  )
})

go_results <- do.call(rbind, results)

# Correct for multiple testing
go_results$adjusted_p_value <- p.adjust(
  go_results$p_value,
  method = "BH"
)

# Sort by adjusted p-value
go_results <- go_results[
  order(go_results$adjusted_p_value),
]

head(go_results)


# Save significant results
significant <- go_results[
  !is.na(go_results$adjusted_p_value) &
    go_results$adjusted_p_value < 0.05,
]

significant

write.table(
  significant,
  "guppy_GO_enrichment_significant.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

head(significant, 20)
