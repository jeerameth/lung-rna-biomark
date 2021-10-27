#--------------------------------------------------------------
# 	link sample id to file id with metadata (BASH)
#--------------------------------------------------------------
awk '{if(/"case_id"/ || /"entity_type": "analyte"/ || /annotation_id/ || (/"file_name"/ && /htseq.counts.gz/)) print $0}' metadata.json > linker1.txt

sed 's/"case_id": "//' linker1.txt >dog.txt
sed 's/"file_name": "//' dog.txt >cat.txt
sed 's/"//g' cat.txt >dog.txt
sed 's/,//g' dog.txt >cat.txt

tr -d " \t" < cat.txt > dog.txt
awk 'NR%2{printf "%s \t",$0;next}{print;}' dog.txt> cat.txt
tr -d " " < cat.txt > linker2.txt

awk '{if(!/anno/) print $0}' linker2.txt > linker3.txt


#--------------------------------------------------------------
# 	link linker to drug
#--------------------------------------------------------------

# Load TCGA  clinical data where 1st and 2nd row were manually deleted

luad <- read.table('clinic/nationwidechildrens.org_clinical_drug_luad.txt', header= TRUE, sep='\t')
lusc <- read.table('clinic/nationwidechildrens.org_clinical_drug_lusc.txt', header= TRUE, sep='\t')

drug <- rbind(luad,lusc)

# Remove unwanted column

drug[28:13] <- NULL
drug[11:7] <- NULL
drug[5:2] <- NULL


# remove n/a and classify to response/noresponse

drug <- drug[!drug$treatment_best_response=="[Unknown]",]
drug <- drug[!drug$treatment_best_response=="[Not Available]",]
drug <- drug[!drug$treatment_best_response=="[Not Applicable]",]

drug$treatment_best_response = gsub("Complete Response","response",drug$treatment_best_response)
drug$treatment_best_response = gsub("Partial Response","response",drug$treatment_best_response)
drug$treatment_best_response = gsub("Stable Disease","response",drug$treatment_best_response)
drug$treatment_best_response = gsub("Clinical Progressive Disease","noresponse",drug$treatment_best_response)


# select drug

drug$pharmaceutical_therapy_drug_name<- tolower(drug$pharmaceutical_therapy_drug_name)
 
# Link case id to file name

linker <- read.table('meta/linker3.txt', header= FALSE, sep='\t')
colnames(linker) <- c('case_id','filename')


# cisplatin !!!!!!!!!!!!!

# cispt <- drug[drug$pharmaceutical_therapy_drug_name=="cisplatin",]
# cisptlink <- merge(linker,cispt,by.x="case_id",by.y="bcr_patient_uuid")
# cisptlink<- cisptlink[!duplicated(cisptlink$case_id),]

# carboplatin !!!!!!!!!!!!!

carbopt <- drug[drug$pharmaceutical_therapy_drug_name=="carboplatin",]
carboptlink <- merge(linker,carbopt,by.x="case_id",by.y="bcr_patient_uuid")
carboptlink<- carboptlink[!duplicated(carboptlink$case_id),]

# check

nrow(cisptlink)
nrow(cisptlink[cisptlink$treatment_best_response=="response",])
nrow(cisptlink[cisptlink$treatment_best_response=="noresponse",])
nrow(carboptlink)
nrow(carboptlink[carboptlink$treatment_best_response=="response",])
nrow(carboptlink[carboptlink$treatment_best_response=="noresponse",])

# https://academic.oup.com/jnci/article/92/3/205/2965042 -> describe responses level (resis = progressive disease)


#--------------------------------------------------------------
# 	DESeq2 + BiomaRt
#--------------------------------------------------------------

# create sample table

directory <- '/home/jeerameth/Desktop/gogobiomarker/htseq'
sampleFiles <- grep("htseq.counts",list.files(path = directory),value=TRUE)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles)

# select drug linker !!!!!!!!!!!!!

# cisptlink$filename <- gsub(".gz","",cisptlink$filename)
# cispttable <- merge(sampleTable ,cisptlink,by.x='fileName',by.y='filename')
# samptable <- cispttable

carboptlink$filename <- gsub(".gz","",carboptlink$filename)
carbopttable <- merge(sampleTable ,carboptlink,by.x='fileName',by.y='filename')
samptable <- carbopttable

# DESeq

library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = samptable,
                                       directory = directory,
                                       design= ~ treatment_best_response)

dds$treatment_best_response
dds$treatment_best_response <- relevel(dds$treatment_best_response, "noresponse")

dds <- DESeq(dds)

resultsNames(dds)
res <- as.data.frame(results(dds))

# BiomaRt (Change Ensembl to HGNC)

library(data.table)
setDT(res, keep.rownames = TRUE)[]
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- res$rn
cutgenes = gsub("\\..*","",genes)
res$cut = cutgenes
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cutgenes,mart= mart)

reshugo <- merge(G_list,res,by.x="ensembl_gene_id",by.y="cut")
reshugo$ensembl_gene_id <- NULL
reshugo$rn <- NULL
reshugo <- reshugo[!reshugo$hgnc_symbol=="",]

write.table(reshugo , file = 'deseq_hugo_drug_respondvnonrespond.txt', quote = FALSE, sep = '\t',row.names=FALSE) 


#--------------------------------------------------------------
# 	GSEA enrichment table (fgsea)
#--------------------------------------------------------------

# retrieve DESeq result  !!!!!!!!!!!!!!!!!!

# reshugo <- read.table('deseq_hugo_cisplatin_respondvnonrespond.txt', header= TRUE, sep='\t')

reshugo <- read.table('deseq_hugo_carboplatin_respondvnonrespond.txt', header= TRUE, sep='\t')

# Gene list prep

library(fgsea)
library(data.table)
library(ggplot2)

geneList = reshugo$log2FoldChange 
names(geneList) = as.character(reshugo$hgnc_symbol)
ranks = sort(geneList, decreasing = TRUE)

# select GO term class !!!!!!!!!!!!!!!!!!!!!!!!!

pathways <- gmtPathways("gotermgmt/c5.go.bp.v7.4.symbols.gmt")

pathways <- gmtPathways("gotermgmt/c5.go.mf.v7.4.symbols.gmt")

pathways <- gmtPathways("gotermgmt/c5.go.cc.v7.4.symbols.gmt")

# GSEA

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf('gseatable_drug_respondvnonrespond_class_top10updown.pdf',width=20, height=5)
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
gseaParam=0.5)
dev.off()

# alternative for top 50 pathway

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=100)
head(fgseaRes)

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=50), pathway]

pdf('gseatable_drug_respondvnonrespond_mf_top50up.pdf',width=20, height=10)
plotGseaTable(pathways[topPathwaysUp], ranks, fgseaRes, 
gseaParam=0.5)
dev.off()

#--------------------------------------------------------------
# 	GSEA (cluster profiler)
#--------------------------------------------------------------

# retrieve DESeq result  !!!!!!!!!!!!!!!!!!

# reshugo <- read.table('deseq_hugo_cisplatin_respondvnonrespond.txt', header= TRUE, sep='\t')

reshugo <- read.table('deseq_hugo_carboplatin_respondvnonrespond.txt', header= TRUE, sep='\t')

# Gene list prep

geneList = reshugo$log2FoldChange 
names(geneList) = as.character(reshugo$hgnc_symbol)
ranks = sort(geneList, decreasing = TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

# select GO class !!!!!!!!!!!!!!!!!!!!!

ego <- gseGO(geneList     = ranks,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              keyType       = 'SYMBOL',
              #nPerm        = 1000, #1000
              minGSSize    = 100, #100
              maxGSSize    = 500, #500
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              seed =FALSE)


pdf('gseadotplot_carboplatin_respondvnonrespond_bp.pdf',width=20, height=10)
dotplot(ego , showCategory=20)
dev.off()


# running score -> select term -> find row number !!!!!!!!!!!!!

ego$Description
library(enrichplot)
match("blood microparticle",ego$Description)

# input rownumber !!!!!!!!!!!!!

pdf('gsearunning_cisplatin_respondvnonrespond_blood_microparticle.pdf',width=10, height=10)
gseaplot2(ego, geneSetID = 2, title = ego$Description[2])
dev.off()

# get leading edge !!!!!!!!!!!!

lead <- ego$core_enrichment[2]
lead <- strsplit(lead, "/")

write.table(lead , file = 'leadingedge_blood_microparticle_cisplatin_respondvnonrespond.txt', quote = FALSE, sep = '\t',row.names=FALSE) 



#--------------------------------------------------------------
# 	ggplot -> violin-like dot plot (usable version is in autobootstrap.r)
#--------------------------------------------------------------



fgseaTable = data.frame(pathway=c('a','b','c','c','d','e','e','e'),NES=c(1,2,4,5,6,6,6,8))



fgseaTable$pathway <- as.factor(fgseaTable$pathway)
head(fgseaTable)

library(ggplot2)
# Basic violin plot
p <- ggplot(fgseaTable, aes(x=NES, y=pathway))
# Set trim argument to FALSE
p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)


pdf('cat.pdf',width=5, height=10)
p
dev.off()


#--------------------------------------------------------------
# 	pre-CARET (gene set selection)
#--------------------------------------------------------------

# Load file !!!!!!!!!!!!!!!!!
df<- read.table('gsea_hallmark_carbopt_table.txt', header= TRUE, sep='\t')
df<- read.table('gsea_hallmark_carbopt_table_ALL.txt', header= TRUE, sep='\t')

# select pathway !!!!!!!!!!!!!!!
pathways <- c('HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
	'HALLMARK_G2M_CHECKPOINT',
	'HALLMARK_E2F_TARGETS')


pathways <- c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
	'HALLMARK_G2M_CHECKPOINT',
	'HALLMARK_E2F_TARGETS')

df <- df[df$pathway %in% pathways,]
lead <- strsplit(df$leadingEdge, split=", ")
names(lead) <- df$pathway


# find common leading edge genes
sumover <- data.frame(gene="a", ratio=0,cat=0)
sumover <- sumover[sumover$gene==0,]

for (i in pathways) {
	sublead <- lead[names(lead)==i]

	total <- vector()
	for (j in 1:length(sublead)){
		total <- append(total,sublead[[j]])
		}
	total

	uniq <- total[!duplicated(total)]
	sum <- data.frame(gene="a", ratio=0)
	sum <- sum[sum$gene==0,]

	for (j in uniq){
	sum <- rbind(sum,data.frame(gene=j, ratio=length(total[total == j])/length(sublead)))
	}
	sum

	sum <- sum[order(-sum$ratio),]

	over <- sum[sum$ratio>0.8,]

	over$cat <- i
	sumover <- rbind(sumover,over)
}
sumover
sumover <- sumover[!duplicated(sumover$gene),]



# create data matrix

directory <- '/home/jeerameth/Desktop/gogobiomarker/htseq'
sampleFiles <- grep("htseq.counts",list.files(path = directory),value=TRUE)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles)

# select drug linker !!!!!!!!!!!!!

# cisptlink$filename <- gsub(".gz","",cisptlink$filename)
# cispttable <- merge(sampleTable ,cisptlink,by.x='fileName',by.y='filename')
# samptable <- cispttable

carboptlink$filename <- gsub(".gz","",carboptlink$filename)
carbopttable <- merge(sampleTable ,carboptlink,by.x='fileName',by.y='filename')
samptable <- carbopttable

# DESeq

library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = samptable,
                                       directory = directory,
                                       design= ~ treatment_best_response)
dds <- DESeq(dds)
count <- counts(dds,normalized=TRUE)
count <- as.data.frame(count)

# BiomaRt

library(data.table)
setDT(count, keep.rownames = TRUE)[]
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- count$rn
cutgenes = gsub("\\..*","",genes)
count$cut = cutgenes
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cutgenes,mart= mart)

counthugo <- merge(G_list,count,by.x="ensembl_gene_id",by.y="cut")
counthugo$ensembl_gene_id <- NULL
counthugo$rn <- NULL
counthugo <- counthugo[!counthugo$hgnc_symbol=="",]

# filter genes
counthugo <- counthugo[counthugo$hgnc_symbol %in% sumover$gene,]


# transpose 
counthugo <- counthugo[!duplicated(counthugo$hgnc_symbol),]

counthugot <- as.data.frame(t(counthugo[,-1]))
colnames(counthugot) <- counthugo$hgnc_symbol 
setDT(counthugot, keep.rownames = TRUE)[]

# select drug linker !!!!!!!!!!!!!
cisptlink$filename <- gsub(".gz","",cisptlink$filename)
cispttable <- merge(counthugot ,cisptlink,by.x='rn',by.y='filename')
countlink <- cispttable

carboptlink$filename <- gsub(".gz","",carboptlink$filename)
carbopttable <- merge(counthugot ,carboptlink,by.x='rn',by.y='filename')
countlink <- carbopttable

# drop column
countlink$case_id <- NULL
countlink$pharmaceutical_therapy_drug_name <- NULL
countlink$rn <- NULL
ncol(countlink) 

countlink$treatment_best_response <- as.factor(countlink$treatment_best_response)

#--------------------------------------------------------------
# 	CARET
#
#	Guide : https://www.machinelearningplus.com/machine-learning/caret-package/
#
#	Flow :
#	1. test subsampling method
#	2. 10x loop for method confirmation
#
#--------------------------------------------------------------



# Data partiton -----------------------------------

library(caret)
str(countlink)
head(countlink[, 1:10])
set.seed(1404)
levels(countlink$treatment_best_response) <- c("noresponse","response")
trainRowNumbers <- createDataPartition(countlink$treatment_best_response, p=0.7, list=FALSE) # try 3070 and 4060
trainData <- countlink[trainRowNumbers,]
testData <- countlink[-trainRowNumbers,]



# CARET #1 (use count of 3 top pathway -> 0.8 common genes) -----------------------------------------------------
set.seed(1404)

imbal_train <- as.data.frame(trainData)
ORIG <- imbal_train
imbal_test  <- as.data.frame(testData)
table(imbal_train$treatment_best_response)

set.seed(1404)
DOWN <- downSample(x = imbal_train[, -ncol(imbal_train)],
                         y = imbal_train$treatment_best_response)
colnames(DOWN)[ncol(DOWN)] <- 'treatment_best_response'
table(DOWN$treatment_best_response)   

set.seed(1404)
UP <- upSample(x = imbal_train[, -ncol(imbal_train)],
                     y = imbal_train$treatment_best_response)
colnames(UP)[ncol(UP)] <- 'treatment_best_response'                         
table(UP$treatment_best_response) 


library(ROSE)

set.seed(1404)
ROSE <- ROSE(treatment_best_response ~ ., data  = imbal_train)$data   
levels(ROSE$treatment_best_response) <- c("noresponse","response")                  
table(ROSE$treatment_best_response) 

#For these data, we’ll use a bagged classification and estimate the area under the ROC curve using five repeats of 10-fold CV.

ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

sampling <- list('ORIG','DOWN','UP')
# sampling <- list('ORIG','DOWN','UP','ROSE')
x= 'DOWN'
test_acc <- function(model) {
  con <- confusionMatrix(predict(model,imbal_test),imbal_test$treatment_best_response)
  con$table
  con$byClass[1:2]
}
test_acct <- function(model) {
  con <- confusionMatrix(predict(model,imbal_train),imbal_train$treatment_best_response)
  con$table
  con$byClass[1:2]
}
allmodels <- data.frame(model="a",sens=0,spec=0)
allmodels <- allmodels[allmodels$model==0,]
#--------------------------------

models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
                    method = "treebag",
                    nbagg = 50,
                    metric = "ROC",
                    trControl = ctrl)
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('treebag_',models_test$rn,sep='') 

allmodels <- rbind(models_test,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))

#-------------------------------- rf

models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
                  method = "earth",
                tuneGrid = data.frame(degree = 2, nprune = 5),
                trControl = ctrl)
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('earth_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))

#-----------------------------------


models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
               method = "rf",
                    metric = "ROC",
                    trControl = ctrl)
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('rf_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))

#-----------------------------------

models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
               method = "svmLinear", 
            trControl = ctrl)
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('svml_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))
#------------------------------------------------------------------------

models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
               method = "glm")
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('glm_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))
#------------------------------------------------------------------------

models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
               method = "rocc")
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('rocc_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))
#------------------------------------------------------------------------

models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
               method = "ranger")
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('ranger_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))
#------------------------------------------------------------------------
models <- list()
fun <- function(x) {
  set.seed(1404)
  model <- train(treatment_best_response ~ ., data = get(x), 
               method = "nnet")
  model <- list(model)
  names(model) <- x
  models <<- c(models,model)
}
lapply(sampling,fun)
models_test <- lapply(models, test_acc)
models_test <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acc), as.vector)))
models_testt <- lapply(models, test_acct)
models_testt <- as.data.frame(do.call("rbind", lapply(lapply(models, test_acct), as.vector)))
colnames(models_testt) <- c("V3","V4")
setDT(models_test, keep.rownames = TRUE)[]
setDT(models_testt, keep.rownames = TRUE)[]
models_test <- merge(models_test,models_testt,by.x="rn",by.y="rn")
models_test$rn <- paste('nnet_',models_test$rn,sep='') 

allmodels <- rbind(allmodels,data.frame(model=models_test$rn,sens_ts=models_test$V1,spec_ts=models_test$V2,sens_tr=models_test$V3,spec_tr=models_test$V4))
allmodels


# 2nd approach using e2071 after realized that linearSVM is bombbbbbb -----------------------------------------


geneList <- c('ISL1','CPS1','COL19A1','LAMB4','CHGA','NECAB2','CTCFL','MAGEA10','SYT5','GDPD2','H19','MYCN','SCN2A','COL2A1','NTRK2','PLAAT5','SYT9','LGALS4','CALB2','OXTR','RPL10P6','ERVH48-1','LERFS','RTL9','LINC00942','KCNJ18')
orimat <- imbal_train[,colnames(imbal_train)%in% geneList]
library("e1071")

X <- as.matrix(orimat[,-ncol(orimat)])
m <- svm(treatment_best_response~., data=orimat, kernel="linear", scale=T,cost='60')
w <- t(m$coefs) %*% X[m$index,]

coef <- t(as.matrix( as.data.frame(coef(m))))
for(column in 1:ncol(sca)){
  coef[, column+1] <- coef[, column+1]/sca[column]
}

p3 <- X %*% coef[,-1] - coef[1]
p3 <- X %*% coef[,-1] -1.651135
p3 <- X %*% coef[,-1]          # <<<<<<< why is it not transposed

coe <- as.matrix(coef[,-1])
p3 <- X %*% coe

cen <- t(as.matrix( as.data.frame(m$x.scale[1])))
# cenmat  <- matrix(, nrow = nrow(imbal_train), ncol = ncol(imbal_train)-1)
# for(column in 1:ncol(cenmat)){
#   cenmat[, column] <- cen[column]
# }

sca <- t(as.matrix( as.data.frame(m$x.scale[2])))

p2 <- predict(m, newdata=imbal_train, decision.values=T)
p1 <- predict(m, newdata=imbal_train, decision.values=F)
p2 <- attr(p2, "decision.values")
# p2 <- X %*% t(w)
p3 <- X %*% t(w) - m$rho
p4 <- X %*% godw - m$rho
p5 <- X 
godp <- data.frame(pred=p1,dec=p2,mod=p3)


# plot heatmap --------------------------------------

df <- table
# rownames(df)<- df[,1]
# df[1] <- NULL
mat <- as.matrix(df)


library(ComplexHeatmap)
library(circlize)


col_col = colorRamp2(seq(min(mat), max(mat), length = 2), c("#EEEEEE", "red"))
font_size =  gpar(fontsize = 10)
pdf('caret_carbotp_roc_sum.pdf',width=4, height=3)
draw(Heatmap(mat, name = "AUC", col = col_col,row_names_gp = font_size ,column_names_gp = font_size,
	, cell_fun = function(j, i, x, y, width, height, fill) 
        {
          grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cluster_columns=FALSE, 
        cluster_rows=FALSE, 

        ))
dev.off()



#PCA ---------------------------------------------------

# create data matrix

directory <- '/home/jeerameth/Desktop/gogobiomarker/htseq'
sampleFiles <- grep("htseq.counts",list.files(path = directory),value=TRUE)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles)

# select drug linker !!!!!!!!!!!!!

cisptlink$filename <- gsub(".gz","",cisptlink$filename)
cispttable <- merge(sampleTable ,cisptlink,by.x='fileName',by.y='filename')
samptable <- cispttable

carboptlink$filename <- gsub(".gz","",carboptlink$filename)
carbopttable <- merge(sampleTable ,carboptlink,by.x='fileName',by.y='filename')
samptable <- carbopttable

# DESeq

library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = samptable,
                                       directory = directory,
                                       design= ~ treatment_best_response)
dds <- DESeq(dds)
count <- counts(dds,normalized=TRUE)
count <- as.data.frame(count)

# BiomaRt

library(data.table)
setDT(count, keep.rownames = TRUE)[]
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- count$rn
cutgenes = gsub("\\..*","",genes)
count$cut = cutgenes
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cutgenes,mart= mart)

counthugo <- merge(G_list,count,by.x="ensembl_gene_id",by.y="cut")
counthugo$ensembl_gene_id <- NULL
counthugo$rn <- NULL
counthugo <- counthugo[!counthugo$hgnc_symbol=="",]

# filter genes
counthugo <- counthugo[!duplicated(counthugo$hgnc_symbol),]

# transpose 

counthugot <- as.data.frame(t(counthugo[,-1]))
colnames(counthugot) <- counthugo$hgnc_symbol 
setDT(counthugot, keep.rownames = TRUE)[]

# select drug linker !!!!!!!!!!!!!
cisptlink$filename <- gsub(".gz","",cisptlink$filename)
cispttable <- merge(counthugot ,cisptlink,by.x='rn',by.y='filename')
countlink <- cispttable

carboptlink$filename <- gsub(".gz","",carboptlink$filename)
carbopttable <- merge(counthugot ,carboptlink,by.x='rn',by.y='filename')
countlink <- carbopttable

# drop column
countlink$case_id <- NULL
countlink$pharmaceutical_therapy_drug_name <- NULL
countlink$rn <- NULL
ncol(countlink) 

countlink$treatment_best_response <- as.factor(countlink$treatment_best_response)


# Data partiton -----------------------------------

library(caret)
str(countlink)
head(countlink[, 1:10])
set.seed(1404)
trainRowNumbers <- createDataPartition(countlink$treatment_best_response, p=0.7, list=FALSE) # try 3070 and 4060
trainData <- countlink[trainRowNumbers,]
testData <- countlink[-trainRowNumbers,]



# CARET #1 (use count of 3 top pathway -> 0.8 common genes) -----------------------------------------------------
set.seed(1404)

imbal_train <- as.data.frame(trainData)
imbal_test  <- as.data.frame(testData)
table(imbal_train$treatment_best_response)


# load filter

test <- countlink[,count]


lungf <- na.omit(lungAll)

lungf <- na.omit(lung)
n <- rownames(lungf)
lungt <- as.data.frame(t(lungf[,-1]))
colnames(lungt) <- n
lungt <- lungt[ , which(apply(lungt, 2, var) != 0)]
lungt$cell_line <- factor(row.names(lungt))
str(lungt) # Check the column types

library(ggfortify)
ncol(lungt)
df <- lungt[1:ncol(lungt)-1]
pca_res <- prcomp(df, scale. = TRUE)

plot <- autoplot(pca_res, data = lungt, colour = 'cell_line')

pdf('PCA.pdf',width=10, height=10)
print(plot)
dev.off()


####################

lungf <- na.omit(lungAll)

write.table(lungf,file ='file_for_PCA.txt',row.names = TRUE,sep = '\t',quote = FALSE)
lungf <- read.delim('file_for_PCA.txt',header=T,sep='\t')


# lungf <- na.omit(lung)
n <- rownames(lungf)
#lungt <- as.data.frame(t(lungf[,-1]))
lungt <- as.data.frame(t(lungf))
colnames(lungt) <- n
lungt <- lungt[ , which(apply(lungt, 2, var) != 0)]
lungt$cell_line <- factor(row.names(lungt))
str(lungt) # Check the column types

library(FactoMineR)
ncol(lungt)
res.pca = PCA(lungt, scale.unit=TRUE, quali.sup=ncol(lungt), graph=F)
plot <- plot.PCA(res.pca, axes=c(1, 2), choix="ind", habillage='ind',label="ind",invisible = "quali")
pdf('PCA_withprim.pdf',width=5, height=5)
print(plot)
dev.off()
var <- as.data.frame(res.pca$var$contrib)
write.table(var,file ='PCA_withprim_dimension_contribution.txt',row.names = TRUE,sep = '\t',quote = FALSE)

