library(data.table);library(ggplot2);library(VennDiagram); library(DESeq2)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=22),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  axis.line = element_line(linewidth = 1.2),axis.ticks=element_line(linewidth=1.2),
  #axis.text.x = element_text(angle =  45,vjust = 0.9, hjust=1)
)
std <- function(x) sd(x,na.rm = T)/sqrt(length(x[!is.na(x)]))
gene_tr_name<-fread("~/Desktop/reference/dmel-all-geneID_transcriptID-r6.04.tr_symbol_sorted",header = F); names(gene_tr_name)<-c("gene_id","gene_name","tr_id","tr_name")


######## Differential expressed gene (DEGs) at each ZT ######## 
hh<-fread("~/Desktop/uorf_ko/rhythm/rhythm_new_nor_dcounts_syq.txt",header = T); nrow(hh) # 9202
zt<-c("ZT0", "ZT4", "ZT8", "ZT12","ZT16", "ZT20")
DE_fem_list <- setNames(vector("list", 6), paste("DE_fem", zt, sep="_"))
DE_mal_list <- setNames(vector("list", 6), paste("DE_mal", zt, sep="_"))
DE_fem_mat<-data.frame(zt0=c(0,0), zt4=c(0,0), zt8=c(0,0), zt12=c(0,0),zt16=c(0,0), zt20=c(0,0)); rownames(DE_fem_mat)<-c("up","down")
DE_mal_mat<-data.frame(zt0=c(0,0), zt4=c(0,0), zt8=c(0,0), zt12=c(0,0),zt16=c(0,0), zt20=c(0,0)); rownames(DE_mal_mat)<-c("up","down")

for (i in 1:length(sel_col)) { #直接运行整个循环可能有与内存读取冲突报错。可以手动赋值 i 运行
  sel_col<-grepl(paste("gene", zt[i], sep = "|"), names(hh))
  hh1<-as.data.frame(hh[, ..sel_col]); rownames(hh1)<-hh1$gene; hh1<-hh1[,-1]
  
  # female
  hh2<-hh1[,grep("_F_",names(hh1))]; hh2[] <- lapply(hh2, as.integer) # 取整数
  colData <- data.frame(line=gsub(".*_(Clk|CS)_.*", "\\1", names(hh2)), 
                        rep=gsub(".*(rep\\d+)$", "\\1", names(hh2))); rownames(colData) <- colnames(hh2)
  hh2_sf <- 1 # The input table has been DESeq2 normalized
  hh2_dds <- DESeqDataSetFromMatrix(countData = hh2, colData, design = ~line) 
  sizeFactors(hh2_dds) <- hh2_sf
  hh2_dds <- estimateDispersions(hh2_dds)
  hh2_dds <- nbinomWaldTest(hh2_dds)
  hh2_res <- as.data.frame(results(hh2_dds, contrast = c('line', 'Clk', 'CS')))
  hh2_res$gene<-rownames(hh2_res); hh2$gene<-rownames(hh2)
  tmp<-merge(hh2_res, hh2, by="gene"); DE_fem_list[[i]]<-tmp; tmp<-as.data.table(tmp)
  DE_fem_mat[1,i]<-nrow(tmp[padj<0.05 & log2FoldChange>= 0.59]); DE_fem_mat[2,i]<-nrow(tmp[padj<0.05 & log2FoldChange<= -0.59])
  
  # male
  hh2<-hh1[,grep("_M_",names(hh1))]; hh2[] <- lapply(hh2, as.integer) # 取整数
  colData <- data.frame(line=gsub(".*_(Clk|CS)_.*", "\\1", names(hh2)), 
                        rep=gsub(".*(rep\\d+)$", "\\1", names(hh2))); rownames(colData) <- colnames(hh2)
  hh2_sf <- 1 # 输入的 table 已经是 DESeq2 normalize好的
  hh2_dds <- DESeqDataSetFromMatrix(countData = hh2, colData, design = ~line) 
  sizeFactors(hh2_dds) <- hh2_sf
  hh2_dds <- estimateDispersions(hh2_dds)
  hh2_dds <- nbinomWaldTest(hh2_dds)
  hh2_res <- as.data.frame(results(hh2_dds, contrast = c('line', 'Clk', 'CS')))
  hh2_res$gene<-rownames(hh2_res); hh2$gene<-rownames(hh2)
  tmp2<-merge(hh2_res, hh2, by="gene"); DE_mal_list[[i]]<-tmp2; tmp2<-as.data.table(tmp2)
  DE_mal_mat[1,i]<-nrow(tmp2[padj<0.05 & log2FoldChange>= 0.59]); DE_mal_mat[2,i]<-nrow(tmp2[padj<0.05 & log2FoldChange<= -0.59])
}

saveRDS(DE_fem_list, "~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/DEGs_fem_eachZT.rds")
saveRDS(DE_mal_list, "~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/DEGs_mal_eachZT.rds")

######## Sleep-related genes in DEGs ######## 
## female Fig.8A
deg_fem<-readRDS("~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/DEGs_fem_eachZT.rds")
sleep<-as.data.frame(fread("~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/sleep_related_535genes.csv"))

for (i in 1:length(deg_fem)){
  x<-as.data.frame(deg_fem[[i]])
  x[is.na(x$padj),]$padj <- 1
  
  x[x$padj < 0.05, ]$gene -> a
  setdiff(x$gene, a) -> b
  intersect(a, sleep[,1]) -> c; setdiff(a, c) -> e
  intersect(b, sleep[,1]) -> d; setdiff(b, d)->f
  print(matrix(c(length(c), length(e), length(d), length(f)), byrow = T, ncol=2))
  cat(i, "all", fisher.test(matrix(c(length(c), length(e), length(d), length(f)), byrow = T, ncol=2))$p.value, "\n")
  print(length(c)/length(e)); print(length(d)/length(f))
}

sleep_ratio_female<-readxl::read_excel("~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/DEG_sleep_fisher.xlsx",sheet = 1,col_names = T)
sleep_ratio_female<-as.data.table(sleep_ratio_female); names(sleep_ratio_female)[1]<-"class"; 
sleep_ratio_female[,c("zt","up_dn"):=.(gsub(".+_zt","",class),gsub("_deg_fem.+","",class))]; 
sleep_ratio_female[,ratio:=Sleep_realted_genes/Non_sleep_realted_genes]
sleep_ratio_female$zt<-factor(sleep_ratio_female$zt,levels = c(0,4,8,12,16,20))
ggplot(sleep_ratio_female[grep("all", class, ignore.case = T)], aes(x = zt, y = ratio, fill=up_dn))+
  geom_bar(width=0.8, position="dodge", stat="identity")+
  # coord_flip()+ #将垂直的条形图变为竖直的条形图，此时注意横坐标的顺序，用 rev 反向一下
  coord_cartesian(ylim = c(0, 0.1)) + 
  labs(x = 'ZT', y = 'Ratio of sleep-related genes', color = NULL)+
  mytheme+theme(legend.position="none",legend.title=element_blank())

## male for Fig.S18A
deg_male<-readRDS("~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/DEGs_mal_eachZT.rds")
sleep<-as.data.frame(fread("~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/sleep_related_535genes.csv"))

for (i in 1:length(deg_fem)){
  x<-as.data.frame(deg_fem[[i]])
  x[is.na(x$padj),]$padj <- 1
  
  x[x$padj < 0.05, ]$gene -> a
  setdiff(x$gene, a) -> b
  intersect(a, sleep[,1]) -> c; setdiff(a, c) -> e
  intersect(b, sleep[,1]) -> d; setdiff(b, d)->f
  print(matrix(c(length(c), length(e), length(d), length(f)), byrow = T, ncol=2))
  cat(i, "all", fisher.test(matrix(c(length(c), length(e), length(d), length(f)), byrow = T, ncol=2))$p.value, "\n")
  print(length(c)/length(e)); print(length(d)/length(f))
}

sleep_ratio_male<-readxl::read_excel("~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/DEG_sleep_fisher.xlsx",sheet = 2,col_names = T)
sleep_ratio_male<-as.data.table(sleep_ratio_male); names(sleep_ratio_male)[1]<-"class"; 
sleep_ratio_male[,c("zt","up_dn"):=.(gsub(".+_zt","",class),gsub("_deg_fem.+","",class))]; 
sleep_ratio_male[,ratio:=Sleep_realted_genes/Non_sleep_realted_genes]
sleep_ratio_male$zt<-factor(sleep_ratio_male$zt,levels = c(0,4,8,12,16,20))
ggplot(sleep_ratio_male[grep("all", class, ignore.case = T)], aes(x = zt, y = ratio, fill=up_dn))+
  geom_bar(width=0.8, position="dodge", stat="identity")+
  # coord_flip()+ #将垂直的条形图变为竖直的条形图，此时注意横坐标的顺序，用 rev 反向一下
  coord_cartesian(ylim = c(0, 0.1)) + 
  labs(x = 'ZT', y = 'Ratio of sleep-related genes', color = NULL)+
  mytheme+theme(legend.position="none",legend.title=element_blank())


######## rhythmiclly expression genes (cycling genes) ######## 
## load FPKM: median FPKM>=1 as expressed gene ####
rpkm<-fread("~/Desktop/uorf_ko/rhythm/rhythm_raw_count_rpkm.csv",header = T)
rpkm1<-melt.data.table(rpkm,id.vars = "gene"); names(rpkm1)<-c("gene","sample","rpkm")
rpkm1[,sample2:=gsub("_rep\\d+","",sample,perl = T)]
rpkm1[,rep_mean:=mean(rpkm),by=.(gene,sample2)] # replicate之间取均值
rpkm2<-unique(rpkm1[,.(gene,sample2,rep_mean)])
rpkm2[,sample3:=gsub("ZT\\d+_","",sample2,perl=T)]; rpkm2<-rpkm2[order(gene, sample3)]
rpkm2[,c("zt_mean","zt_median"):=.(mean(rep_mean),median(rep_mean)),by=.(gene,sample3)] # 每个品系内，同一基因在不同时间点之间取均值/中位数
rpkm2[,c("zt_max","zt_min"):=.(max(rep_mean),min(rep_mean)),by=.(gene,sample3)]
rpkm2[,fc:= (zt_max+0.01)/(zt_min+0.01)] # amplitud=max/min
rpkm2<-rpkm2[order(gene,sample3)]

rpkm1_fwt<-unique(rpkm2[sample3=="CS_F" & zt_median>=1, .(gene, sample3, zt_mean, zt_median, zt_max, zt_min, fc)]) # 8788
rpkm1_fmut<-unique(rpkm2[sample3=="Clk_F" & zt_median>=1, .(gene, sample3, zt_mean, zt_median, zt_max, zt_min, fc)]) # 8857
rpkm1_mwt<-unique(rpkm2[sample3=="CS_M" & zt_median>=1, .(gene, sample3, zt_mean, zt_median, zt_max, zt_min, fc)]) # 8920
rpkm1_mmut<-unique(rpkm2[sample3=="Clk_M" & zt_median>=1, .(gene, sample3, zt_mean, zt_median, zt_max, zt_min, fc)]) # 8992

## load MetaCycle results ####
fwt<-fread("~/Desktop/uorf_ko/rhythm/RNA-Seq/MetaCycle/rhythm_new_nor_dcounts_MetaCyc_WT_fem_syq.csv",header = T)
fmut<-fread("~/Desktop/uorf_ko/rhythm/RNA-Seq/MetaCycle/rhythm_new_nor_dcounts_MetaCyc_Mut_fem_syq.csv",header = T)
mwt<-fread("~/Desktop/uorf_ko/rhythm/RNA-Seq/MetaCycle/rhythm_new_nor_dcounts_MetaCyc_WT_male_syq.csv",header = T)
mmut<-fread("~/Desktop/uorf_ko/rhythm/RNA-Seq/MetaCycle/rhythm_new_nor_dcounts_MetaCyc_Mut_male_syq.csv",header = T)

## both qvalues <= 0.05 & median RPKM>=1 & max/min fc >= 1.5 ####
# female WT vs female KO
fwt_005_15<- fwt[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_fwt[fc>=1.5]$gene]; nrow(fwt_005_15) # 440
fmut_005_15<- fmut[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_fmut[fc>=1.5]$gene]; nrow(fmut_005_15) # 648
venn.plot <-venn.diagram(list(unique(fwt_005_15$gene),unique(fmut_005_15$gene)),category.names=c("WT_F","Mut_F"),resolution=600,filename=NULL,lwd=1.5,cex=4,cat.pos=c(0,0),cat.fontface=2,na="remove",cat.cex=2,fill=rainbow(2))
grid.newpage();grid.draw(venn.plot) # wt, common, mut: 151, 289, 359
# male WT vs male KO
mwt_005_15<- mwt[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_mwt[fc>=1.5]$gene]; nrow(mwt_005_15) # 472
mmut_005_15<- mmut[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_mmut[fc>=1.5]$gene]; nrow(mmut_005_15) # 512
venn.plot <-venn.diagram(list(unique(mwt_005_15$gene),unique(mmut_005_15$gene)),category.names=c("WT_M","Mut_M"),resolution=600,filename=NULL,lwd=1.5,cex=4,cat.pos=c(0,0),cat.fontface=2,na="remove",cat.cex=2,fill=rainbow(2))
grid.newpage();grid.draw(venn.plot) # wt, common, mut: 242, 230, 282

rhythm_p005_fc15<-list(fwt=fwt_005_15, fmut=fmut_005_15, mwt=mwt_005_15, mmut=mmut_005_15)
saveRDS(rhythm_p005_fc15, "~/Desktop/uorf_ko/rhythm/manuascript/20240909 plos Biology/rhythm_p005_fc15.rds")

##  Overlap between female and male ####
## diffrent qvalue cutoff: cycle gene 和 sex difference ####
# both qvalues <= 0.05 & median RPKM>=1 & max/min fc >= 1.5
# female WT vs female KO
fwt_005_15<- fwt[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_fwt[fc>=1.5]$gene]; nrow(fwt_005_15) # 440
fmut_005_15<- fmut[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_fmut[fc>=1.5]$gene]; nrow(fmut_005_15) # 648
venn.plot <-venn.diagram(list(unique(fwt_005_15$gene),unique(fmut_005_15$gene)),category.names=c("WT_F","Mut_F"),resolution=600,filename=NULL,lwd=1.5,cex=4,cat.pos=c(0,0),cat.fontface=2,na="remove",cat.cex=2,fill=rainbow(2))
grid.newpage();grid.draw(venn.plot) # wt, common, mut: 151, 289, 359
# male WT vs male KO
mwt_005_15<- mwt[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_mwt[fc>=1.5]$gene]; nrow(mwt_005_15) # 472
mmut_005_15<- mmut[ARS_BH.Q<=0.05 & JTK_BH.Q<=0.05 & gene%in%rpkm1_mmut[fc>=1.5]$gene]; nrow(mmut_005_15) # 512
venn.plot <-venn.diagram(list(unique(mwt_005_15$gene),unique(mmut_005_15$gene)),category.names=c("WT_M","Mut_M"),resolution=600,filename=NULL,lwd=1.5,cex=4,cat.pos=c(0,0),cat.fontface=2,na="remove",cat.cex=2,fill=rainbow(2))
grid.newpage();grid.draw(venn.plot) # wt, common, mut: 242, 230, 282

# both qvalues <= 0.01 & median RPKM>=1 & max/min fc >= 1.5 
# female WT vs female KO
fwt_001_15<- fwt[ARS_BH.Q<=0.01 & JTK_BH.Q<=0.01 & gene%in%rpkm1_fwt[fc>=1.5]$gene]; nrow(fwt_001_15) # 306
fmut_001_15<- fmut[ARS_BH.Q<=0.01 & JTK_BH.Q<=0.01 & gene%in%rpkm1_fmut[fc>=1.5]$gene]; nrow(fmut_001_15) # 463
venn.plot <-venn.diagram(list(unique(fwt_001_15$gene),unique(fmut_001_15$gene)),category.names=c("WT_F","Mut_F"),resolution=600,filename=NULL,lwd=1.5,cex=4,cat.pos=c(0,0),cat.fontface=2,na="remove",cat.cex=2,fill=rainbow(2))
grid.newpage();grid.draw(venn.plot) 
# male WT vs male KO
mwt_001_15<- mwt[ARS_BH.Q<=0.01 & JTK_BH.Q<=0.01 & gene%in%rpkm1_mwt[fc>=1.5]$gene]; nrow(mwt_001_15) # 292
mmut_001_15<- mmut[ARS_BH.Q<=0.01 & JTK_BH.Q<=0.01 & gene%in%rpkm1_mmut[fc>=1.5]$gene]; nrow(mmut_001_15) # 293
venn.plot <-venn.diagram(list(unique(mwt_001_15$gene),unique(mmut_001_15$gene)),category.names=c("WT_M","Mut_M"),resolution=600,filename=NULL,lwd=1.5,cex=4,cat.pos=c(0,0),cat.fontface=2,na="remove",cat.cex=2,fill=rainbow(2))
grid.newpage();grid.draw(venn.plot)

## barplot showing number of cycling genes under different qvalue ####
# female
qdata<-data.table(num=c(440,648,306,463))
qdata[,qvalue:=c("q<=0.05","q<=0.05","q<=0.01","q<=0.01")]; qdata$qvalue<-factor(qdata$qvalue,levels = c("q<=0.05","q<=0.01"))
qdata[,geno:=rep(c("WT","KO"),times=2)]; qdata$geno<-factor(qdata$geno,levels = c("WT","KO"))
ggplot(qdata, aes(x = qvalue, y = num, fill=factor(geno)))+
  geom_bar(width=0.8,position="dodge", stat="identity",color="black",linewidth=1.1)+
  scale_fill_manual(values = c("gray", "#0F99B2"))+
  labs(x = 'Significance thresholds ', y = 'Counts', color = NULL)+
  mytheme+theme(legend.position="none",legend.title=element_blank())

# male
qdata<-data.table(num=c(472,512,292,293))
qdata[,qvalue:=c("q<=0.05","q<=0.05","q<=0.01","q<=0.01")]; qdata$qvalue<-factor(qdata$qvalue,levels = c("q<=0.05","q<=0.01"))
qdata[,geno:=rep(c("WT","KO"),times=2)]; qdata$geno<-factor(qdata$geno,levels = c("WT","KO"))
ggplot(qdata, aes(x = qvalue, y = num, fill=factor(geno)))+
  geom_bar(width=0.8,position="dodge", stat="identity",color="black",linewidth=1.1)+
  scale_fill_manual(values = c("gray", "#0F99B2"))+
  labs(x = 'Significance thresholds ', y = 'Counts', color = NULL)+
  mytheme+theme(legend.position="none",legend.title=element_blank())

save.image("~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/Clock_uORF_code/Fig8.RData")

