library(plyr)
library(stringr)
library(GEOquery)

###### Create phenotypic data set
gse <- getGEO("GSE52529")
e1 <- gse$`GSE52529-GPL11154_series_matrix.txt.gz`
e2 <- gse$`GSE52529-GPL16791_series_matrix.txt.gz`
pd <- rbind.fill(pData(e1),pData(e2))

# extract SRA number
pd$SRA_Experiment <- unlist(lapply(str_split(as.character(pd$relation.1),"=",n=2), function(x) x[2]))

# load srp
srp <- read.csv("extdata/SraRunInfo_SRP033135.csv")
srpsmall <- srp[,c("SampleName", "Run", "Experiment", "Sample", "BioSample",  
                   "avgLength", "SampleType", "LibraryLayout", 
                   "Model", "ScientificName", "LoadDate", "ReleaseDate")]
colnames(srpsmall)[which(colnames(srpsmall) == "Experiment")] <- "SRA_Experiment"
coldata <- merge(pd, srpsmall, by ="SRA_Experiment", all.x = TRUE)
coldata$download_path <- paste0(file.path(coldata$supplementary_file_1, coldata$Run, coldata$Run), ".sra")
rownames(coldata) <- coldata$Run

# add batch effects
firstline <- read.table("summaryStats/fastqBatches.txt")
fqInfo <- array(NA, dim = c(length(firstline$V2), 7))
dat1 <- data.frame(ldply(str_split(firstline$V2[!grepl("C268PACXX", firstline$V2)], ":", n=7)))
colnames(dat1) <- c("instrument", "runID", "fcID", "fcLane", "tile","xtile", "ytile")
dat2 <- data.frame(ldply(str_split(firstline$V2[grepl("C268PACXX", firstline$V2)], ":", n=5))) 
colnames(dat2) <- c("fcID", "fcLane", "tile","xtile", "ytile")
fqDat1 <- cbind(data.frame("Run" = str_sub(firstline$V1[!grepl("C268PACXX", firstline$V2)], start = 2, end = -5)), dat1)
fqDat2 <- cbind(data.frame("Run" = str_sub(firstline$V1[grepl("C268PACXX", firstline$V2)], start = 2, end = -5)), dat2)
fqDat <- rbind.fill(fqDat1, fqDat2)
coldata <- merge(coldata, fqDat, by = "Run")
rownames(coldata) <- coldata$Run

# subset for relevant phenotypic data
coldata$hour <- coldata$characteristics_ch1.1
levels(coldata$hour) <- list("hour72" = "hour post serum-switch: 72", 
                             "hour0" = "hour post serum-switch: 0", 
                             "Hour24" ="hour post serum-switch: 24", 
                             "Hour48" = "hour post serum-switch: 48")

coldata$sampleType <- factor(ifelse(grepl("Bulk", coldata$characteristics_ch1.5), "bulk",
                                    "SC"), levels=c("SC", "bulk"))
coldata$control <- factor(coldata$characteristics_ch1.3)
coldata$debris <- factor(coldata$characteristics_ch1.2)
coldata$numcells <- factor(coldata$characteristics_ch1.4)

coldata$sampleName <- str_sub(coldata$title, start = 6)
coldata <- coldata[, colnames(coldata) %in% 
                  c("Run", "geo_accession", "Model", "instrument", "runID", "fcID", 
                    "fcLane", "hour", "sampleType", "sampleName", "source_name_ch1",
                    "control", "debris", "numcells", "description")]


####### Create ExpressionSet object
esetNorm <- read.delim("fpkms/GSE52529_fpkm_matrix.txt", header = TRUE, sep = "\t")
esetNormBulk <- read.delim("fpkms/GSE52529_truseq_fpkm_matrix.txt", header = TRUE, sep = "\t")
colnames(esetNormBulk) <- unlist(lapply(str_split(colnames(esetNormBulk), "_"), 
              function(x){ paste(x[[1]], as.numeric(x[[2]]) + 1, sep="_")}))
eset <- cbind(esetNorm, esetNormBulk)
eset <- eset[, match(as.character(coldata$description), colnames(eset))] # order
colnames(eset) <- coldata$Run

# Create ExpressionSet
trapnell2014myoblasthuman = ExpressionSet(assayData = as.matrix(eset), 
                             phenoData = AnnotatedDataFrame(coldata))

# Save ExpressionSet
save(trapnell2014myoblasthuman, file="data/trapnell2014myoblasthuman.rda")
