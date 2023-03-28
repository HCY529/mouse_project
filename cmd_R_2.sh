library(DSS)
data.1 = read.table("SRR833527_1_val_1_bismark_bt2_pe.bismark.cov.gz",header = T)
data.2 = read.table("SRR833528_1_val_1_bismark_bt2_pe.bismark.cov.gz",header = T)
BSobj <- makeBSseqData(list(data.1,data.2),c("C1","N1"))[1:1000,]
BSobj =  makeBSseqData(list(data.1,data.2),c("C1","N1"))[1:1000,]
data.1 = read.table("SRR833527_1_val_1_bismark_bt2_pe.bismark.cov.gz",header = T,check.names=F)
data.2 = read.table("SRR833528_1_val_1_bismark_bt2_pe.bismark.cov.gz",header = T,check.names=F)
BSobj =  makeBSseqData(list(data.1,data.2),c("C1","N1"))[1:1000,]
data.1 = read.table("SRR833527_1_val_1_bismark_bt2_pe.bismark.cov",header = T)
data.2 = read.table("SRR833528_1_val_1_bismark_bt2_pe.bismark.cov",header = T)
BSobj =  makeBSseqData(list(data.1,data.2),c("C1","N1"))[1:1000,]
data.1
data.1 = read.table("SRR833527_1_val_1_bismark_bt2_pe.bismark.cov",header = F)
data.1
data.2 = read.table("SRR833528_1_val_1_bismark_bt2_pe.bismark.cov",header = F)
BSobj =  makeBSseqData(list(data.1,data.2),c("C1","N1"))[1:1000,]
data.2 = read.table("SRR833528_1_val_1_bismark_bt2_pe.dss.input.txt",header = T)
data.1 = read.table("SRR833527_1_val_1_bismark_bt2_pe.dss.input.txt",header = T)
BSobj =  makeBSseqData(list(data.1,data.2),c("C1","N1"))[1:1000,]
dmlTest <- DMLtest(BSobj, group1=c("C1"), group2=c("N1"))
dmlTest <- DMLtest(BSobj, group1=c("C1"), group2=c("N1"),smoothing=TRUE)
dmls <- callDML(dmlTest.sm, p.threshold=0.001)
head(dmlTest.sm)
dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
showOneDMR(dmrs[1,], BSobj)
write.table(dmrs,file=paste0("dmrs.txt"),sep = "\t",append = FALSE,row.names = FALSE, col.names = TRUE, quote = FALSE)
library(ChIPseeker)
BiocManager::install("CHIPseeker")
library(ChIPseeker)
BiocManager::install("ChIPseeker")
library(ChIPseeker)
BiocManager::install("ChIPseeker")
BiocManager::install("GenomicFeatures")
BiocManager::install("curl")
BiocManager::install("GO.db")
BiocManager::install("GenomicFeatures",force = TRUE)
keggplot <- function(i,gene_list,showCategory,outfile,p){
    kk <- enrichKEGG(gene_list[[i]],organism='dosa',keyType = 'kegg',                
                 pvalueCutoff=p, pAdjustMethod='none',                 
                 qvalueCutoff=p)
    print(barplot(kk,showCategory = showCategory,title = ""))
    print(dotplot(kk,showCategory = showCategory,title = ""))
    print(heatplot(kk,showCategory = showCategory))
    write.table(data.frame(kk@result),file = outfile,sep="\t", row.names = FALSE, col.names =TRUE, quote =FALSE) #output 
}
goplot <- function(i,gene_list,showCategory,outfile,p){
    ego.BP <- enrichGO(gene_list[[i]],OrgDb = rice, keyType = "RAP", ont="BP",              
                 pvalueCutoff=p, pAdjustMethod='none',                 
                 qvalueCutoff=p)
    ego.CC <- enrichGO(gene_list[[i]],OrgDb = rice, keyType = "RAP", ont="CC",              
                 pvalueCutoff=p, pAdjustMethod='none',                 
                 qvalueCutoff=p)
    ego.MF <- enrichGO(gene_list[[i]],OrgDb = rice, keyType = "RAP", ont="MF",              
                 pvalueCutoff=p, pAdjustMethod='none',                 
                 qvalueCutoff=p)
    ego <- enrichGO(gene_list[[i]],OrgDb = rice, keyType = "RAP", ont="ALL",              
                 pvalueCutoff=p, pAdjustMethod='none',                 
                 qvalueCutoff=p)
    write.table(data.frame(ego@result),file = outfile,sep="\t", row.names = FALSE, col.names =TRUE, quote =FALSE) #output 
    print(barplot(ego.BP,showCategory = showCategory,title = paste0("BP")))
    print(dotplot(ego.BP,showCategory = showCategory,title = paste0("BP")))
    print(heatplot(ego.BP,showCategory = showCategory))
    print(plotGOgraph(ego.BP))
    print(barplot(ego.CC,showCategory = showCategory,title = paste0("CC")))
    print(dotplot(ego.CC,showCategory = showCategory,title = paste0("CC")))
    print(heatplot(ego.CC,showCategory = showCategory))
    print(plotGOgraph(ego.CC))
    print(barplot(ego.MF,showCategory = showCategory,title = paste0("MF")))
    print(dotplot(ego.MF,showCategory = showCategory,title = paste0("MF")))
    print(heatplot(ego.MF,showCategory = showCategory))
    print(plotGOgraph(ego.MF))
}
lapply(1, keggplot,showCategory=10,gene_list=transcript,outfile=paste0("dmrs_kegg.txt"),p=1)
library(clusterProfiler)
BiocManager::install("clusterProfiler")
q()
