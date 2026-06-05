
#install.packages("BiocManager", lib='~/bin/R/library',repos='https://cloud.r-project.org')
#library(BiocManager, lib='~/bin/R/library')
#BiocManager::install(version = "3.18", lib='~/bin/R/library')
#BiocManager::install("biomaRt", lib = "~/bin/R/library")
library('biomaRt', lib='~/bin/R/library')
library('topGO', lib='~/bin/R/library')

aln_block_length=1000
aln_qual=40
#biomart_attributes=c("ensembl_gene_id","entrezgene_id","uniprot_gn_id","kegg_enzyme","chromosome_name","start_position","end_position")
#active_mart=guppy
biomart_attributes=c("ensembl_gene_id","external_gene_name", "entrezgene_id","uniprot_gn_id","go_id","name_1006","namespace_1003", "chromosome_name","start_position","end_position")

guppy <- useEnsembl(
  biomart = "genes",
  dataset = "preticulata_gene_ensembl"
)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1] # First argument is the PAF file

#aln<-read.table("candidate_window_sel1.paf",fill=T)
aln <- read.table(input_file, fill=TRUE)

aln<-aln[aln$V11 > aln_block_length & aln$V12 >= aln_qual,]

regions<-paste0(aln$V6,":",as.integer(aln$V8),":",as.integer(aln$V9))

regions2 <- gsub("^chr", "LG", regions)

tmp_biomart<-getBM(attributes = biomart_attributes,
                   filters= "chromosomal_region",
                   values=regions2,
                   mart=guppy)


output_file <- paste0(input_file, ".out")

#write.table(tmp_biomart, "biomart_sel1.txt", sep = "\t", quote=F)
write.table(tmp_biomart, output_file, sep = "\t", quote=F)

#Get all guppy gene IDs (universe)
universe_biomart <- getBM(
  attributes = c("ensembl_gene_id", "go_id"),
  mart = guppy
)
geneUniverse <- unique(universe_biomart$ensembl_gene_id)

#Get interesting genes from above
int_genes <- tmp_biomart$ensembl_gene_id

# 3. Prepare gene list for topGO
geneList <- factor(as.integer(geneUniverse %in% int_genes))
names(geneList) <- geneUniverse

#Prepare gene2GO mapping
gene2GO <- tapply(universe_biomart$go_id, universe_biomart$ensembl_gene_id, function(x) unique(na.omit(x)))

#Make topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x) x == 1,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)

#Fisher's exact test for enrichment
resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

#Summarize results
allRes <- GenTable(GOdata, classicFisher = resultFisher)
output_go <- paste0(input_file, "GO.out")

#write.table(allRes, "GO_enrichment_results.txt", sep="\t", row.names=F, quote=F)
write.table(allRes, output_go, sep = "\t", row.names=F, quote=F)

print(allRes)
