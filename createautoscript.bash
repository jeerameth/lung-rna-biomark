#!/bin/bash

# sed -i -e 's/\r$//' createautoscript.bash

for drug in 'carbopt'
do

	if [ $drug = 'carbopt' ]
		then
			echo 'if(file.exists("deseq_'$drug'_table_10set.txt")) {
			deseqTable <- read.table("deseq_'$drug'_table_10set.txt", header= TRUE, sep="\t")
			}else{' > autobootstrap.r

		else
			echo 'if(file.exists("deseq_'$drug'_table_10set.txt")) {
			deseqTable <- read.table("deseq_'$drug'_table_10set.txt", header= TRUE, sep="\t")
			}else{' >> autobootstrap.r
		fi

	echo 'luad <- read.table("clinic/nationwidechildrens.org_clinical_drug_luad.txt", header= TRUE, sep="\t")
	lusc <- read.table("clinic/nationwidechildrens.org_clinical_drug_lusc.txt", header= TRUE, sep="\t")
	drug <- rbind(luad,lusc)
	drug[28:13] <- NULL
	drug[11:7] <- NULL
	drug[5:2] <- NULL
	drug <- drug[!drug$treatment_best_response=="[Unknown]",]
	drug <- drug[!drug$treatment_best_response=="[Not Available]",]
	drug <- drug[!drug$treatment_best_response=="[Not Applicable]",]
	drug$treatment_best_response = gsub("Complete Response","response",drug$treatment_best_response)
	drug$treatment_best_response = gsub("Partial Response","response",drug$treatment_best_response)
	drug$treatment_best_response = gsub("Stable Disease","response",drug$treatment_best_response)
	drug$treatment_best_response = gsub("Clinical Progressive Disease","noresponse",drug$treatment_best_response)
	drug$pharmaceutical_therapy_drug_name<- tolower(drug$pharmaceutical_therapy_drug_name)
	linker <- read.table("meta/linker3.txt", header= FALSE, sep="\t")
	colnames(linker) <- c("case_id","filename")
	cispt <- drug[drug$pharmaceutical_therapy_drug_name=="cisplatin",]
	cisptlink <- merge(linker,cispt,by.x="case_id",by.y="bcr_patient_uuid")
	cisptlink<- cisptlink[!duplicated(cisptlink$case_id),]
	carbopt <- drug[drug$pharmaceutical_therapy_drug_name=="carboplatin",]
	carboptlink <- merge(linker,carbopt,by.x="case_id",by.y="bcr_patient_uuid")
	carboptlink<- carboptlink[!duplicated(carboptlink$case_id),]
	nrow(cisptlink)
	nrow(cisptlink[cisptlink$treatment_best_response=="response",])
	nrow(cisptlink[cisptlink$treatment_best_response=="noresponse",])
	nrow(carboptlink)
	nrow(carboptlink[carboptlink$treatment_best_response=="response",])
	nrow(carboptlink[carboptlink$treatment_best_response=="noresponse",])
	directory <- "/home/jeerameth/Desktop/gogobiomarker/htseq"
	sampleFiles <- grep("htseq.counts",list.files(path = directory),value=TRUE)
	sampleTable <- data.frame(sampleName = sampleFiles,
	                          fileName = sampleFiles)' >> autobootstrap.r


	echo $drug'link$filename <- gsub(".gz","",'$drug'link$filename)
	'$drug'table <- merge(sampleTable ,'$drug'link,by.x="fileName",by.y="filename")
	samptable <- '$drug'table' >> autobootstrap.r


	echo $drug'link$filename <- gsub(".gz","",'$drug'link$filename)
	'$drug'table <- merge(sampleTable ,'$drug'link,by.x="fileName",by.y="filename")
	samptable <- '$drug'table
	samptable <- samptable[order(samptable$treatment_best_response),]
	row.names(samptable) <- NULL
	}' >> autobootstrap.r


	for i in {1..10}
	do
		echo 'if(!file.exists("deseq_'$drug'_table_10set.txt")) {' >> autobootstrap.r
		
		if [ $drug = 'cispt' ]
		then
			echo 'nresp <- sample(1:7, 5, replace=F)
		resp <- sample(8:54, 5, replace=F)
		ransamptable <- samptable[c(nresp,resp),]' >> autobootstrap.r

		else
			echo 'nresp <- sample(1:15, 15, replace=F)
		resp <- sample(16:47, 15, replace=F)
		ransamptable <- samptable[c(nresp,resp),]' >> autobootstrap.r
		fi

		echo 'library("DESeq2")
		dds <- DESeqDataSetFromHTSeqCount(sampleTable = ransamptable,
		                                       directory = directory,
		                                       design= ~ treatment_best_response)
		dds$treatment_best_response <- relevel(dds$treatment_best_response, "noresponse")
		dds <- DESeq(dds)
		res <- as.data.frame(results(dds))



		library(data.table)
		setDT(res, keep.rownames = TRUE)[]
		library("biomaRt")
		mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
		genes <- res$rn
		cutgenes = gsub("\\..*","",genes)
		res$cut = cutgenes
		G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cutgenes,mart= mart)

		reshugo <- merge(G_list,res,by.x="ensembl_gene_id",by.y="cut")
		reshugo$ensembl_gene_id <- NULL
		reshugo$rn <- NULL
		reshugo <- reshugo[!reshugo$hgnc_symbol=="",]
		geneList = reshugo$log2FoldChange 
		names(geneList) = as.character(reshugo$hgnc_symbol)
		reshugoDEG = reshugo[(reshugo$log2FoldChange>1 |reshugo$log2FoldChange< -1) & reshugo$pval<0.05,]
		geneListDEG = reshugoDEG$log2FoldChange 
		names(geneListDEG) = as.character(reshugoDEG$hgnc_symbol)


		}else{
			geneList <- deseqTable[(deseqTable$rn == "log2FoldChange" & deseqTable$iteration == '$i'),2:ncol(deseqTable)]
			geneList$iteration <- NULL
			name <- colnames(geneList)
			geneList <- as.numeric(geneList)
			names(geneList) = as.character(name)


			deseqTableOne <- deseqTable[(deseqTable$iteration == '$i'),]
			deseqTableOne$iteration <- NULL
			reshugo <- as.data.frame(t(deseqTableOne[,-1]))
			colnames(reshugo) <- deseqTableOne$rn
			reshugoDEG = reshugo[(reshugo$log2FoldChange > 1 |reshugo$log2FoldChange < -1) & reshugo$pval < 0.05,]
			library(data.table)
			setDT(reshugoDEG, keep.rownames = TRUE)[]
			reshugoDEG <- as.data.frame(reshugoDEG)
			reshugoDEG$rn = gsub("\\..*","",reshugoDEG$rn)
			reshugoDEG <- reshugoDEG[(reshugoDEG$rn != "NA"),]
			reshugoDEG <- as.data.frame(reshugoDEG)
			geneListDEG = reshugoDEG$log2FoldChange 
			names(geneListDEG) = as.character(reshugoDEG$rn)
		}


		library(fgsea)
		library(data.table)
		library(ggplot2)
		
		ranks = sort(geneList, decreasing = TRUE)
		pathways <- gmtPathways("gotermgmt/h.all.v7.4.symbols.gmt")
		fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
		head(fgseaRes)
		fgseaRes <- as.data.frame(fgseaRes)
		fgseaCut <- fgseaRes[fgseaRes$padj<=0.05,]
		fgsea'$i' <- fgseaCut$pathway
		fgseaCut[2:5] <- NULL


		ranksDEG = sort(geneListDEG, decreasing = TRUE)
		pathways <- gmtPathways("gotermgmt/h.all.v7.4.symbols.gmt")
		fgseaResDEG <- fgsea(pathways, ranksDEG, minSize=15, maxSize=500)
		head(fgseaResDEG)
		fgseaResDEG <- as.data.frame(fgseaResDEG)
		fgseaCutDEG <- fgseaResDEG[fgseaResDEG$padj<=0.05,]
		fgseaDEG'$i' <- fgseaCutDEG$pathway
		fgseaCutDEG[2:5] <- NULL


		library(clusterProfiler)
		library(org.Hs.eg.db)

		ego <- enrichGO(gene          = names(geneListDEG),
		                universe      = names(geneList),
		                keyType       = "SYMBOL",
		                OrgDb         = org.Hs.eg.db,
		                ont           = "ALL",
		                pAdjustMethod = "fdr",
		                pvalueCutoff  = 0.05,)
		enrich <- summary(ego)
		enrichCut <- enrich[enrich$pvalue<=0.05,]
		enrichCut$term <- paste("GO_" ,paste(enrichCut$ONTOLOGY ,paste(":" ,enrichCut$Description,sep=""),sep=""),sep="")
		enrichCut$neglogpval <- -log10(enrichCut$pval)
		enrichCut[1:10] <- NULL



if(!file.exists("deseq_'$drug'_table_10set.txt")) {

		reshugo <- reshugo[!duplicated(reshugo$hgnc_symbol),]
		reshugot <- as.data.frame(t(reshugo[,-1]))
		colnames(reshugot) <- reshugo$hgnc_symbol 
		setDT(reshugot, keep.rownames = TRUE)[]
		reshugot$iteration <- '$i >> autobootstrap.r

		if [ $i == 1 ]
		then
			echo 'deseqTable <-  as.data.frame(reshugot)
			deseqTable' >> autobootstrap.r

		else
			echo 'library("dplyr")
			deseqTable <- bind_rows(deseqTable, as.data.frame(reshugot))' >> autobootstrap.r
		fi

		echo '}' >> autobootstrap.r


		if [ $i == 1 ]
		then
			echo 'fgseaTable <- fgseaCut
			fgseaTableDEG <- fgseaCutDEG
			enrichTable <- enrichCut' >> autobootstrap.r

		else
			echo 'fgseaTable <- rbind(fgseaTable,fgseaCut)
			fgseaTableDEG <- rbind(fgseaTableDEG ,fgseaCutDEG )
			enrichTable <- rbind(enrichTable ,enrichCut )' >> autobootstrap.r
		fi

		echo ' cat("***********************************

				Finish iteration '$i' of '$drug'.

				*************************************")' >> autobootstrap.r

	done




	echo 'fgseaTable$pathway <- as.factor(fgseaTable$pathway)
	head(fgseaTable)

	library(ggplot2)
	cat("****************************************plot1")
	p <- ggplot(fgseaTable, aes(x=NES, y=pathway))
	p <- p + geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5,stackratio=0.2) + labs(title="Enriched GSs in '$drug' sensitive lung cancer")
	pdf("gsea_hallmark_'$drug'_dotplot_stacked.pdf",width=10, height=10)
	p
	dev.off()

	library(ggplot2)
	cat("****************************************plot2")
	p <- ggplot(fgseaTableDEG, aes(x=NES, y=pathway))
	p <- p + geom_dotplot(binaxis="y", stackdir="center", dotsize=0.5,stackratio=0.2) + labs(title="Enriched GSs in '$drug' sensitive lung cancer")
	pdf("gsea_hallmark_'$drug'_DEG_dotplot_stacked.pdf",width=10, height=10)
	p
	dev.off()

	uniqTable <- fgseaTable[!duplicated(fgseaTable$pathway),]
	sumTable <- data.frame(pathway="a", mNES=0,ratio=0)
	sumTable <- sumTable[sumTable$pathway==0,]

	for (i in uniqTable$pathway){
	sumTable <- rbind(sumTable,data.frame(pathway=i, mNES=mean(fgseaTable[fgseaTable$pathway == i,2]),ratio=nrow(fgseaTable[fgseaTable$pathway == i,])/10))
	}
	sumTable
	sumTable <- sumTable[sumTable$ration>=0.7,]


	uniqTableDEG <- fgseaTableDEG[!duplicated(fgseaTableDEG$pathway),]
	sumTableDEG <- data.frame(pathway="a", mNES=0,ratio=0)
	sumTableDEG <- sumTableDEG[sumTableDEG$pathway==0,]

	for (i in uniqTableDEG$pathway){
	sumTableDEG <- rbind(sumTableDEG,data.frame(pathway=i, mNES=mean(fgseaTableDEG[fgseaTableDEG$pathway == i,2]),ratio=nrow(fgseaTableDEG[fgseaTableDEG$pathway == i,])/10))
	}
	sumTableDEG

	cat("****************************************plot3")
	p <- ggplot(sumTable, aes(x=mNES, y=pathway, size=ratio))
	p <- p + geom_point()
	pdf("gsea_hallmark_'$drug'_dotplot_summed.pdf",width=8, height=6)
	p
	dev.off()

	sumTable$type <- "Transcriptome"
	sumTableDEG$type <- "DEG"
	allTable <- rbind(sumTable,sumTableDEG)

	cat("****************************************plot4")

	p <- ggplot(allTable, aes(x=mNES, y=pathway, size=ratio, color=type))
	p <- p + geom_point()
	pdf("gsea_hallmark_'$drug'_dotplot_summed_withDEG.pdf",width=8, height=6)
	p
	dev.off()




	uniqTableEn <- enrichTable[!duplicated(enrichTable$term),]
	sumTableEn <- data.frame(term="a", neglogpval=0)
	sumTableEn <- sumTableEn[sumTableEn$term==0,]
	for (i in uniqTableEn$term){
	sumTableEn <- rbind(sumTableEn,data.frame(term=i, neglogpval=mean(enrichTable[enrichTable$term == i,2]),ratio=nrow(enrichTable[enrichTable$term == i,])/10))
	}
	sumTableEn <- sumTableEn[sumTableEn$ratio>=0.8,]
	sumTableEn


	cat("****************************************plot5")

	library(ggplot2)
	p <- ggplot(sumTableEn, aes(x=neglogpval, y=term, size=ratio))
	p <- p + geom_point() + labs(title="Over represented GO terms in '$drug' sensitive lung cancer")
	pdf("goterm_'$drug'_DEG_dotplot_stacked.pdf",width=20, height=10)
	p
	dev.off()
		cat("**************************************** finished plot")'  >> autobootstrap.r






	echo 'fgseaTable$leadingEdge <- vapply(fgseaTable$leadingEdge, paste, collapse = ", ", character(1L))
	cat("**************************************** 1 chunk")
	write.table(fgseaTable, file = "gsea_hallmark_'$drug'_table.txt", quote = FALSE, sep = "\t",row.names=FALSE)
	cat("**************************************** 2 chunk")
	if(!file.exists("deseq_'$drug'_table_10set.txt")) {
	write.table(deseqTable, file = "deseq_'$drug'_table_10set.txt", quote = FALSE, sep = "\t",row.names=FALSE) 
	}
	cat("**************************************** 3 chunk")' >> autobootstrap.r







	echo 'df <- data.frame(vec1 = rep(NA, max(sapply(list(fgsea1' >> autobootstrap.r


	for i in {2..10}
	do
		echo ',fgsea'$i >> autobootstrap.r

	done

	echo '),length))))' >> autobootstrap.r


	for i in {1..10}
	do
		echo 'df[1:length(fgsea'$i'), '$i'] <- fgsea'$i >> autobootstrap.r

	done
	
	echo 'write.table(df, file = "gsea_hallmark_'$drug'_10set.txt", quote = FALSE, sep = "\t",row.names=FALSE)
	cat("**************************************** last chunk") ' >> autobootstrap.r



done
