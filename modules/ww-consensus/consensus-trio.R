#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is one argument: if not, return an error
if (length(args) != 4) {
  stop("4 arguments must be provided, GATK file, SAMtools file, Mutect file and the base file name.n", call.=FALSE)
}

GATKfile <- args[1]
SAMfile <- args[2]
MutectFile <- args[3]
baseName <- args[4]

## Load Libraries
library(tidyverse)
library(dplyr)

print("Process GATK variants")
G <- read.delim(file = GATKfile, row.names = NULL,
                  header=FALSE, stringsAsFactors = FALSE, skip = 1)
annHead <- read.delim(file = GATKfile, row.names = NULL,
                        nrows = 1, header=FALSE, stringsAsFactors = FALSE)
extraHead <- c("QUAL1", "DP.GATK", "CHR", "POS", "END", "REF",
               "ALT", "QUAL.GATK", "FILTER", "INFO.GATK", "FORMAT.GATK", "MAGIC.GATK")
colnames(G) <- c(annHead[1,], extraHead)
G[, "QUAL1"]<- NULL
G$DP.GATK <- as.integer(G$DP.GATK)
G$VariantID <- paste(G$CHR, G$Gene.refGene, G$POS, G$REF, G$ALT, sep = "-")
G$AD.GATK <- as.integer(gsub("^[^:]+:[^,]+,|:.*", "", G$MAGIC.GATK));
G$VAF.GATK <- 100*as.numeric(G$AD.GATK)/as.numeric(G$DP.GATK);
G <- Filter(function(y)!all(y == "."), G);
G[G == ""] <- NA
G <- G %>% dplyr::mutate_if(is.factor, as.character);
G <- G %>% mutate_if(is.integer, as.character);


print("Process SAMTOOLS variants")
S <- read.delim(file = SAMfile, row.names=NULL, 
                  header = FALSE, skip = 1)
annHead <- read.delim(file = SAMfile, row.names = NULL, 
                        nrows = 1, header=FALSE, stringsAsFactors = FALSE)
extraHead <- c("QUAL1", "DP.SAM", "CHR", "POS", "END", "REF", 
               "ALT", "QUAL.SAM", "FILTER", "INFO.SAM", "FORMAT.SAM", "MAGIC.SAM")
colnames(S) <- c(annHead[1,], extraHead)
#DP4 = Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles,
# used in variant calling. Sum can be smaller than DP because low-quality bases are not counted.
S[,"QUAL1"]<- NULL
S$VariantID <- paste(S$CHR, S$Gene.refGene, S$POS, S$REF, S$ALT, sep = "-")
a <- gsub(".*:", "",S$MAGIC.SAM);
REF <- as.numeric(gsub(",.*$","",a));
S$AD.SAM <- as.numeric(gsub("^.*,","",a));
S$DP.SAM <- REF+S$AD.SAM;
S$VAF.SAM <- 100*S$AD.SAM/S$DP.SAM;
S <- Filter(function(y)!all(y == "."), S);
S[S == ""] <- NA
S <- S %>% mutate_if(is.factor, as.character);
S <- S %>% mutate_if(is.integer, as.character);


print("Process Mutect variants")

Mu <- read.delim(file = MutectFile, row.names=NULL, 
                header = FALSE, skip = 1)
annHead <- read.delim(file = MutectFile, row.names = NULL, 
                      nrows = 1, header=FALSE, stringsAsFactors = FALSE)
extraHead <- c("col1", "DP.Mu", "CHR", "POS", "END", "REF", 
               "ALT", "col2", "FILTER.Mu", "INFO.Mu", "FORMAT.Mu", "MAGIC.Mu")
colnames(Mu) <- c(annHead[1,], extraHead)
Mu[,"col1"]<- NULL
Mu$VariantID <- paste(Mu$CHR, Mu$Gene.refGene, Mu$POS, Mu$REF, Mu$ALT, sep = "-")
Mu$DP.Mu <- as.integer(Mu$DP.Mu)
Mu$VariantID <- paste(Mu$CHR, Mu$Gene.refGene, Mu$POS, Mu$REF, Mu$ALT, sep = "-")
Mu$AD.Mu <- as.integer(gsub("^[^:]+:[^,]+,|:.*", "", Mu$MAGIC.Mu));
Mu$VAF.Mu <- 100*as.numeric(Mu$AD.Mu)/as.numeric(Mu$DP.Mu);
Mu <- Filter(function(y)!all(y == "."), Mu);
Mu[Mu == ""] <- NA
Mu <- Mu %>% dplyr::mutate_if(is.factor, as.character);
Mu <- Mu %>% mutate_if(is.integer, as.character);

## Deal with flubbed annotation issues where the same variant is now annotated two different ways!  ARGH!
commonCols <- c("VariantID","CHR", "POS", "REF", "ALT", "Chr", "Start", "End", "Ref", "Alt")
GATKCols <- colnames(G)[!colnames(G) %in% c(colnames(S), colnames(Mu))]
SAMCols <- colnames(S)[!colnames(S) %in% c(colnames(G), colnames(Mu))]
MuCols <- colnames(Mu)[!colnames(Mu) %in% c(colnames(G), colnames(S))]

GMerge <- G %>% select(c(commonCols, GATKCols))
SMerge <- S %>% select(c(commonCols, SAMCols))
MuMerge <- Mu %>% select(c(commonCols, MuCols))

variants <- full_join(full_join(GMerge, SMerge), MuMerge)
variants$Type <- ifelse(nchar(variants$REF) == nchar(variants$ALT), "SNV", "INDEL");

reannotate <- left_join(variants, G, by = c("VariantID", "CHR", "POS", "REF", "ALT"))
reannotate <- left_join(reannotate, S, by = c("VariantID", "CHR", "POS", "REF", "ALT"))
reannotate <- left_join(reannotate, Mu, by = c("VariantID", "CHR", "POS", "REF", "ALT"))

print("DEBUG: Columns in reannotate after joins:")
print(colnames(reannotate))
print(paste("DEBUG: nrow =", nrow(reannotate)))
print("DEBUG: AD columns check:")
print(paste("AD.GATK exists:", "AD.GATK" %in% colnames(reannotate)))
print(paste("AD.SAM exists:", "AD.SAM" %in% colnames(reannotate)))
print(paste("AD.Mu exists:", "AD.Mu" %in% colnames(reannotate)))
if ("AD.GATK" %in% colnames(reannotate)) {
  print("DEBUG: AD.GATK values:")
  print(reannotate$AD.GATK)
}

reannotate$Confidence <- case_when(
  !is.na(reannotate$AD.GATK) & !is.na(reannotate$AD.SAM) & !is.na(reannotate$AD.Mu) ~ "conftier1",
  (is.na(reannotate$AD.GATK) | is.na(reannotate$AD.Mu)) & !is.na(reannotate$AD.SAM) & reannotate$Type == "SNV" ~ "conftier2",
  !(is.na(reannotate$AD.GATK) | is.na(reannotate$AD.Mu)) & is.na(reannotate$AD.SAM) & reannotate$Type == "INDEL" ~ "conftier2",
  TRUE ~ "conftier3"
)
write.table(reannotate, file = paste0(baseName, ".consensus.tsv"), sep = "\t",
            row.names = F)

