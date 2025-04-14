library(data.table);library(ggplot2);library(VennDiagram)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=22),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  axis.line = element_line(linewidth = 1.2),axis.ticks=element_line(linewidth=1.2),
  #axis.text.x = element_text(angle =  45,vjust = 0.9, hjust=1)
)
#### BLS caluculation ####
library(ape)
tree_fly <- read.tree(file="~/Desktop/uorf_ko/uORF_gain_loss/dm6-droGri2_iNodes.nh")

jj_bls<-fread("~/Desktop/uorf_ko/uORF_gain_loss/uORF_matrix_dm6_PacBioSim_27species_update_noCDS.csv",header = T)
names(jj_bls)[1]<-"id";jj_bls[,c(-4,-26:-29)]->jj_bls # remove old sim and non-drosophila species
jj_bls[,presence:=rowSums(jj_bls[,2:24])]; jj_bls1<-jj_bls[presence>0,] # remove uorfs not present in any of 23 flies
jj_bls1[,pos_id:=paste(seqname,gene_pos,sep="_")] # genomic position
jj_bls1[,id_mel:=paste(transcriptID,position_withoutgap+1,sep="_")] # transcript position in dmel

names(jj_bls1)[3]<-"droSim1" #为了不改 nh 文件里的物种名字，这里改成 droSim1与 nh 文件一致，后续在改回来
list_fly <- apply(jj_bls1[,2:24],1,function(x){
  sp <- which(x>0)
  paste0(colnames(jj_bls1)[2:24][sp],collapse=",")
})
list_fly <- as.data.frame(list_fly)
colnames(list_fly) <- c("sp")
list_fly[,1] <- as.character(list_fly[,1])

# a long time for BLS calculation
list_fly$score <- 0
for(i in 1:nrow(list_fly)){
  list_fly$score[i] <- sum(keep.tip(tree_fly,unlist(strsplit(list_fly[i,1],",")))$edge.length)
}

total_score<-sum(tree_fly$edge.length) #6.557502
list_fly$BLS <- list_fly$score/total_score
names(list_fly)[1]<-"present_sp"
list_fly$present_sp<-gsub("droSim1","PacBioSim",list_fly$present_sp)
cbind(jj_bls1,list_fly)->jj_bls2
fwrite(jj_bls2,"~/Desktop/uorf_ko/uORF_gain_loss/uORF_matrix_dm6_PacBioSim_23Droso_update_noCDS.BLS.txt",col.names = T,sep="\t",quote = F)

#### violin plot of drosophila clock uATGs distribution for Fig. 1A and S2A ####
pro_gene<-fread("~/Desktop/reference/dmel-all-geneID_transcriptID-r6.04.protein_symbol_sorted",header = F)
pro_gene<-as.data.table(unique(pro_gene$V1)); names(pro_gene)<-"gene_id"; nrow(pro_gene) # 13917 protein-coding gene_ID
clock<-fread("~/Desktop/uorf_ko/rhythm/FlyBase_CircadianRhythmIDs.txt",header = F) 
names(clock)<-"gene_id"; nrow(clock) # 169 clock genes
core_clock<-readxl::read_excel("~/Desktop/uorf_ko/rhythm/clock_uorf_statistics.xlsx",sheet = 2)
names(core_clock)<-c("gene_name","gene_id"); nrow(core_clock) # 10

hh1<-fread("~/Desktop/uorf_ko/uORF_gain_loss/uORF_matrix_dm6_PacBioSim_23Droso_update_noCDS.BLS.txt",header = T) # 自己算的，步骤在这一部分uORF BLS based on 23 drosophila species
# hh1<-hh[!is.na(geneID)&dm6==1,]; nrow(hh1) # some tr were not annotated in currently used GTF; total 70713 redundant uATGs present in dmel
# hh1[,gene_pos_id:=paste(geneID,seqname,gene_pos,sep="_")] 
# hh1<-hh1[,-2:-24]

## all tr considered, but only use genomic-unique uATGs
hh2<-hh1[!duplicated(gene_pos_id)]; nrow(hh2); nrow(hh1[!duplicated(pos_id)]) # gene_genomic unique vs  genomic unique: 37311 vs 36590
hh2[,uorf_num:=.N,by=.(geneID)]
hh3<-unique(hh2[,.(geneID,uorf_num)]) # 7607个基因至少有一个 uorf
# 要把不含有 uorf 的基因也考虑进来，记为 0。而不是只比较有 uorf 的基因
gene_uorf<-merge(pro_gene,hh3,by.x="gene_id",by.y="geneID",all.x=T)
gene_uorf[is.na(gene_uorf)]<-0
gene_uorf[,class:=ifelse(gene_id%in%clock$gene_id,"Clock","Others")]
table(gene_uorf$class); sum(gene_uorf[class=="Clock"]$uorf_num); sum(gene_uorf[class=="Others"]$uorf_num)
# clock: 152, 1137; others: 13756, 36174
fisher.test(matrix(c(152,1137,13756,36174),nrow = 2,byrow =T))$`p.value` # 5.78007e-42

summary(gene_uorf[class=="Clock"]$uorf_num); summary(gene_uorf[class=="Others"]$uorf_num)
wilcox.test(gene_uorf[class=="Clock"]$uorf_num,gene_uorf[class=="Others"]$uorf_num,alternative = "greater")$`p.value` # 2.638656e-23

gene_uorf$class<-factor(gene_uorf$class,levels = c("Clock","Others"))
## most tr considered
ggplot(gene_uorf,aes(x=class,y=uorf_num))+
  geom_violin(aes(fill=class, color=class),width=1.3,adjust = 2.5, linewidth=1.2)+
  geom_boxplot(width=0.05,outlier.size = 1,fatten = 3, lwd=1.2)+
  # scale_color_manual(name = NULL,values = c(RColorBrewer::brewer.pal(3, 'Dark2')[2], 'grey40'),breaks = c('Clock', 'Others'),labels = c('Circadian genes, 152', 'Other genes, 13,756')) +
  scale_color_manual(values = c("#E04643", "#686F76"))+
  scale_fill_manual(values = alpha(c("#E04643", "#686F76"), 0.8))+
  labs(x = ' ', y = 'uORF number') + 
  coord_cartesian(ylim = c(0, 30), expand = 0) +
  mytheme + theme(legend.position="none")

tmp<-copy(gene_uorf); names(tmp)<-c("Gene_id","uORF_count", "Gene_categories")
fwrite(tmp,"~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/fig1A_data.csv", sep=",", col.names = T,row.names = F,quote = F)

## most tr considered: Fig. S2A
most_tr<-fread("~/Desktop/uorf_ko/rhythm/most_tr_all_smp_kallisto_tpm.mergedRep.txt",header = T)
jj1<-hh1[transcriptID%in%most_tr$tr]; nrow(jj1) # 19696个uATG; 先选 tr 再去重;
jj2<-jj1[!duplicated(gene_pos_id)]; nrow(jj2) # 19696 uATG; 因为每个基因只选了一个 tr，所以没有 uATG 重叠
jj2[,uorf_num:=.N,by=.(geneID)]
jj3<-unique(jj2[,.(geneID,uorf_num)]) # 6440 uATG; 要把不含有 uorf 的基因也考虑进来，而不是只比较有 uorf 的基因
gene_uorf<-merge(pro_gene,jj3,by.x="gene_id",by.y="geneID",all.x=T)
gene_uorf[is.na(gene_uorf)]<-0
gene_uorf[,class:=ifelse(gene_id%in%clock$gene_id,"Clock","Others")]
table(gene_uorf$class)
sum(gene_uorf[class=="Clock"]$uorf_num); sum(gene_uorf[class=="Others"]$uorf_num)
# Clock: 152, 505; others:13765, 19191
fisher.test(matrix(c(152,505,13765,19191),nrow = 2,byrow =T))$`p.value` # 3.443007e-23
summary(gene_uorf[class=="Clock"]$uorf_num); summary(gene_uorf[class=="Others"]$uorf_num)
wilcox.test(gene_uorf[class=="Clock"]$uorf_num,gene_uorf[class=="Others"]$uorf_num,alternative = "greater")$`p.value` # 1.443477e-14

gene_uorf$class<-factor(gene_uorf$class,levels = c("Clock","Others"))
ggplot(gene_uorf,aes(x=class,y=uorf_num))+
  geom_violin(aes(fill=class, color=class),width=1.3,adjust = 2.5, linewidth=1.2)+
  geom_boxplot(width=0.05,outlier.size = 1,fatten = 3, lwd=1.2)+
  # scale_color_manual(name = NULL,values = c(RColorBrewer::brewer.pal(3, 'Dark2')[2], 'grey40'),breaks = c('Clock', 'Others'),labels = c('Circadian genes, 152', 'Other genes, 13,756')) +
  scale_color_manual(values = c("#E04643", "#686F76"))+
  scale_fill_manual(values = alpha(c("#E04643", "#686F76"), 0.8))+
  labs(x = ' ', y = 'uORF number') + 
  coord_cartesian(ylim = c(0, 20), expand = 0) +
  mytheme + theme(legend.position="none")

tmp<-copy(gene_uorf); names(tmp)<-c("Gene_id","uORF_count", "Gene_categories")
fwrite(tmp,"~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/figS2A_data.csv", sep=",", col.names = T,row.names = F,quote = F)

#### Proportion of highly conserved uORFs in  different gene categories for Fig. 1B ####
cc<-hh1[!duplicated(pos_id)] # total 36590 genomic-unique uATGs in dmel
cc[,gene_class:= ifelse(geneID%in%core_clock$gene_id, "core", ifelse(geneID%in%clock$gene_id, "nonCore_clock","others"))]
table(cc$gene_class)  # total 82, 1045, 35428 genomic-unique uATGs in core-clock, nonCore-clock and other genes
table(cc[BLS==1]$gene_class) # total 7, 22, 359 highly conserved genomic-unique uATGs in core-clock, nonCore-clock and other genes
#### The translational efficiency (TE) of different gene classes for Fig. 1C ####
te_fh<-fread("~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/fig1c_data_FH.csv",header = T)
ggplot(te_fh, aes(x = Gene_class, y = TE))+
  geom_jitter(data = subset(te_fh[!is.infinite(TE)], Gene_class != "Others"),
              aes(color = Gene_class), size=2.2, alpha = 0.7, position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = Gene_class, group=factor(Gene_class)), position = position_dodge(width = 0.85), outlier.shape = NA, fill = '#00000000',
               width = 0.7, lwd=1.5 , fatten = 1.5)+ 
  scale_color_manual(values = c('#FC8767', '#108EA9', 'black', 'gray'))+
  coord_cartesian(ylim = c(0, 4)) + labs(x = '', y = 'Translation efficiency (TE)') + 
  mytheme + theme(legend.position="none") 

te_mh<-fread("~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/fig1c_data_MH.csv",header = T)
ggplot(te_mh, aes(x = Gene_class, y = TE))+
  geom_jitter(data = subset(te_mh[!is.infinite(TE)], Gene_class != "Others"),
              aes(color = Gene_class), size=2.2, alpha = 0.7, position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = Gene_class, group=factor(Gene_class)), position = position_dodge(width = 0.85), outlier.shape = NA, fill = '#00000000',
               width = 0.7, lwd=1.5 , fatten = 1.5)+ 
  scale_color_manual(values = c('#FC8767', '#108EA9', 'black', 'gray'))+
  coord_cartesian(ylim = c(0, 4)) + labs(x = '', y = 'Translation efficiency (TE)') + 
  mytheme + theme(legend.position="none") 

#### BLS for uATGs in different gene categories for Fig.S1B and S2B ####
figS1B<-fread("~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/figS1B_data.csv", header = T)
ggplot(figS1B, aes(x = BLS, color = Gene_category)) + 
  stat_ecdf(linewidth = 1.2) +
  scale_color_manual(name = NULL,
                     values = c(RColorBrewer::brewer.pal(3, 'Dark2')[1:2], 'grey40'),
                     breaks = c('Core_clock', 'Clock', 'Others'),
                     labels = c('Core clock genes (7)', 'Non-core clock genes (120)', 'Other genes (7,480)')) +
  labs(x = 'BLS', y = 'ECDF') + 
  coord_cartesian(xlim = c(0, 1), expand = 0) +
  mytheme + theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = c(0.99, 0.01), legend.justification = c(1, 0))


figS2B<-fread("~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/figS2B.csv", header = T)
ggplot(figS2B, aes(x = BLS, color = Gene_category)) +
  stat_ecdf(linewidth = 1.2) +
  scale_color_manual(name = NULL,
                     values = c(RColorBrewer::brewer.pal(3, 'Dark2')[1:2], 'grey40'),
                     breaks = c('Core_clock', 'Clock', 'Others'),
                     labels = c('Core clock genes (10)', 'Non-core clock genes (135)', 'Other genes (11,511)')) +
  labs(x = 'BLS', y = 'ECDF') + 
  coord_cartesian(xlim = c(0, 1), expand = 0) +
  mytheme + theme(plot.margin = unit(c(1,1,1,1), "cm"), legend.position = c(0.99, 0.01), legend.justification = c(1, 0))

save.image("~/Desktop/uorf_ko/rhythm/manuascript/20250403 plos biology/Clock_uORF_code/Fig1.RData")




