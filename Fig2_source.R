setwd("~/syq/project/uorf/gbrowser_graph/")
###################################################################################################
## loading package
###################################################################################################
suppressPackageStartupMessages(library('GenomicFeatures'))
suppressPackageStartupMessages(library('rtracklayer'))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
library(RColorBrewer)
library(readxl)


###################################################################################################
## make transcriptome Txdb and uORF annotation
###################################################################################################
fly6.04 <- loadDb('./dmel-all-r6.04.TxDb')
fly.exon <- exonsBy(fly6.04, by = 'tx', use.names = TRUE)
fly.genes <- genes(fly6.04)

# # prepare unique complete uORF annotations
# uorf.bed12 <- fread('../results/all_aug_orfs_complete.bed')
# uorf.bed12.uniq <- unique(
#     uorf.bed12, by = c('chrom', 'start', 'end', 'strand', 'blocksize', 'blockstart'))
# fwrite(uorf.bed12.uniq, 'all_aug_orfs_complete_uniq.bed', sep = '\t', quote = F)
# # tail -n +2 all_aug_orfs_complete_uniq.bed \
# #     | bedtools bed12tobed6 -i stdin >all_aug_orfs_complete_uniq_6col.bed

uorf.annot <- fread('all_aug_orfs_complete_6col.bed') #保留所有的 uORF，不用像 ZH 上面做的那样 unique
setnames(uorf.annot, names(uorf.annot),
         c('chrm', 'start', 'end', 'name', 'transcript', 'strand'))

fly.isoform <- as.data.frame(transcriptsBy(fly6.04))
setDT(fly.isoform)
setnames(fly.isoform, c('group', 'group_name'), c('gene_id', 'gene_name'))

fly.cds <- cdsBy(fly6.04, 'tx', use.names = TRUE)
fly.utr5 <- fiveUTRsByTranscript(fly6.04, use.names = TRUE)
fly.utr3 <- threeUTRsByTranscript(fly6.04, use.names = TRUE)

###################################################################################################
## prepare bw files and helper funcitons for visualizing RNA-seq coverage
###################################################################################################
bw.files <- scan('BigWigFileList_syq.txt', what = '')
bw.files.meta <- cbind(
  matrix(bw.files[1:26], ncol = 2, byrow = TRUE),
  matrix(bw.files[27:52], ncol = 2, byrow = TRUE)
)

rownames(bw.files.meta) <- c(
  'mature_oocyte', 'em_0_2h', 'em_2_6h', 'em_6_12h', 'em_12_24h', 'larva', 'pupa',
  'female_body1', 'female_body2', 'male_body1', 'male_body2', 'female_head', 'male_head')
colnames(bw.files.meta) <- c(
  'mrna.forward', 'mrna.reverse', 'ribo.forward', 'ribo.reverse')


ExtractGeneCov <- function(gene, samples = rownames(bw.files.meta)){
  # get gene range
  g.intv <- fly.genes[gene]
  chrm <- as.character(runValue(seqnames(g.intv)))
  g.strand <- as.character(runValue(strand(g.intv)))
  g.intv <- keepSeqlevels(g.intv, chrm)
  # gt bw file list
  if(g.strand == '+'){
    bw.files.list <- bw.files.meta[samples, c(1, 3)]
  }else{
    bw.files.list <- bw.files.meta[samples, c(2, 4)]
  }
  bw.files.list <- BigWigFileList(as.vector(bw.files.list))
  # extract coverage
  mat.cov <- sapply(bw.files.list, function(bw){
    bw <- import(bw, which = g.intv)
    bw.cov <- coverage(keepSeqlevels(bw, chrm), weight = 'score')
    as.vector(Views(bw.cov[[1]], ranges(g.intv))[[1]])
  })
  if(g.strand == '-'){
    mat.cov <- mat.cov[nrow(mat.cov):1, ] * (-1)
  }
  mat.cov
}


ExtractGeneInfo <- function(gene, genome.ind = FALSE, drop.isoform = NULL){
  g.exons.df <- tryCatch(
    {
      g.isoform <- fly.isoform[gene_name == gene, tx_name]
      as.data.frame(fly.exon[names(fly.exon) %in% g.isoform])[2:7]
    },
    error = function(e){warning('Gene not found!'); NULL}
  )
  if(is.null(g.exons.df)) return(g.exons.df)
  setnames(g.exons.df, c('group_name', 'seqnames'), c('tx_name', 'chrm'))
  
  if(!genome.ind){
    scols <- c('start', 'end')
    g.strand <- as.character(g.exons.df$strand[1])
    if(g.strand == '+'){
      m <- min(g.exons.df$start, g.exons.df$end)
      g.exons.df[scols] <- g.exons.df[scols] - m + 1
    }else if(g.strand == '-'){
      m <- max(g.exons.df$start, g.exons.df$end)
      g.exons.df[scols] <- m - g.exons.df[rev(scols)] + 1
    }
  }
  if(!is.null(drop.isoform))
    g.exons.df <- g.exons.df[!g.exons.df$tx_name %in% drop.isoform, ]
  
  g.exons.df
}


RemoveIntron <- function(df.cov, df.exons, keep.ori = FALSE){
  df.cov <- as.data.table(df.cov)
  df.exons <- as.data.table(df.exons)[, .(start, end)]
  res <- df.cov[posn %inrange% df.exons]
  res <- as.data.frame(res)
  if(!keep.ori){
    res$posn <- seq_len(nrow(res))
  }
  res
}


MakeCovPlot <- function(g, ..., plt.samples = rownames(bw.files.meta),
                        log.scale = FALSE, nm = FALSE, plot = TRUE){
  #g <- 'FBgn0034876'
  g.sample <- c(t(outer(c('mrna', 'ribo'), plt.samples, FUN = paste, sep = '|')))
  g.mat <-ExtractGeneCov(g, plt.samples)
  if(log.scale){
    g.mat <- apply(g.mat, 2, function(x){
      m <- min(x[x > 0])
      log(x/m + 1)
    })
  }
  colnames(g.mat) <- g.sample
  
  g.df <- as.data.frame(g.mat)
  g.df$posn <- seq_len(nrow(g.df))
  
  g.exons <- ExtractGeneInfo(g, ...)
  g.plt <- RemoveIntron(g.df, g.exons)
  
  g.plt <- melt(
    g.plt, id.vars = 'posn', variable.name = 'track', value.name = 'coverage')
  if(nm){
    setDT(g.plt)
    g.plt[, coverage := coverage / max(coverage), by = .(track)]
    setDF(g.plt)
  }
  if(is.data.table(g.plt)){
    print(head(g.plt))
  }
  g.plt[c('type', 'sample')] <- tstrsplit(g.plt$track, '|', fixed = TRUE)
  # could also use levels = plt.samples
  g.plt$sample <- factor(g.plt$sample, levels = plt.samples)
  
  if(plot){
    ggplot(g.plt, aes(x = posn, y = coverage)) +
      facet_grid(sample ~ type, scales = 'free_y') +
      geom_ribbon(aes(ymax = coverage, ymin = 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_classic()
  }else{
    invisible(g.plt)
  }
}

###################################################################################################
## plot gene model with geom_ribbon
###################################################################################################
## helper functions
PlotGeneModel <- function(gene){
  g.intv <- fly.genes[gene]
  ggbio::autoplot(fly6.04, wh = range(g.intv))
}


MakeFeatTable <- function(gene){
  g.intv <- fly.genes[gene]
  chrm <- as.character(seqnames(g.intv)[1])
  g.strand <- as.character(strand(g.intv)[1])
  g.intv <- keepSeqlevels(g.intv, chrm)
  
  g.isoform <- fly.isoform[gene_name == gene]
  
  # cds/utr intervals
  df.utr5 <- as.data.frame(fly.utr5[names(fly.utr5) %in% g.isoform$tx_name])
  df.cds <- as.data.frame(fly.cds[names(fly.cds) %in% g.isoform$tx_name])
  df.utr3 <- as.data.frame(fly.utr3[names(fly.utr3) %in% g.isoform$tx_name])
  scols <- c('group_name', 'seqnames', 'start', 'end', 'strand')
  df <- list(utr5 = df.utr5[scols], cds = df.cds[scols], utr3 = df.utr3[scols])
  df <- rbindlist(df, idcol = 'type')
  setnames(df, 'seqnames', 'chrm')
  
  df.uorf <- uorf.annot[
    transcript %in% g.isoform$tx_name,
    .(type = 'uorf', group_name = name, chrm, start, end, strand)]
  rbind(df, df.uorf)
}

# change log:
#  2017-12-19: add `xlim` and `drop.manual` to customize plots
GeneModelExon <- function(gene, drop.isoform = NULL, xlim = NULL, drop.manual = NULL){
  gene.feat <- MakeFeatTable(gene)
  gene.exons.raw <- ExtractGeneInfo(gene, genome.ind = TRUE)
  
  # convert genome indices into gene indices
  scols <- c('start', 'end')
  if(gene.exons.raw$strand[1] == '+'){
    start.ind <- min(gene.exons.raw$start)
    #set(gene.feat, NULL, scols, gene.feat[, scols, with = F] - start.ind + 1)
    gene.feat[,c("start","end"):=.(start-start.ind + 1,end-start.ind + 1)]
  }else if(gene.exons.raw$strand[1] == '-'){
    start.ind <- max(gene.exons.raw$end)
    #set(gene.feat, NULL, scols, start.ind - gene.feat[, rev(scols), with = F] + 1)
    gene.feat[,c("start","end"):=.(start.ind - end + 1,start.ind - start + 1)]
    
  }
  # map current range boundary to exon only indices,
  # and drop isoforms redundant isforms if nececcary
  gene.drop <- drop.isoform
  if(!is.null(drop.isoform))
    gene.exons <- ExtractGeneInfo(gene, drop.isoform = gene.drop)
  else
    gene.exons <- ExtractGeneInfo(gene)
  m <- inrange(
    seq_len(max(gene.exons[scols])),
    lower = gene.exons$start, upper = gene.exons$end)
  m <- cumsum(m)
  gene.feat <- gene.feat[!group_name %in% gene.drop]
  gene.feat[, `:=`(start = m[start], end = m[end])]
  
  # only plot gene model in the specified window
  if(!is.null(xlim)){
    gene.feat <- gene.feat[start <= xlim[2] & end >= xlim[1] ]
    gene.feat[, end := pmin(xlim[2], end)]
    gene.feat[, start := pmax(xlim[1], start)]
  }
  # manual drop some redundant features
  
  if(!is.null(drop.manual)){
    gene.feat <- gene.feat[!group_name %in% drop.manual]
  }
  gene.feat
}


GeneModelExonPlot <- function(gene, ...){
  # extract gene feature
  df <- GeneModelExon(gene, ...)
  
  # group feature by feature type: transcript (cds, utr), then uORF
  df <- df[order(type == 'uorf', group_name, decreasing = TRUE)]
  df[, y := rleid(group_name)]
  # UTR and CDS have different height
  df[, ymin := ifelse(type == 'cds', y - 0.4, y - 0.25)]
  df[, ymax := ifelse(type == 'cds', y + 0.4, y + 0.25)]
  # each row is a unique group
  df[, fid := factor(1:.N)]
  # reshape start:end to make a rectangle
  plt <- melt(df, measure.vars = c('start', 'end'), value.name = 'posn')
  
  # do not discriminate UTR types in color mode
  plt[type %in% c('utr3', 'utr5'), type := 'utr']
  
  ggplot(plt, aes(x = posn)) +
    # concatenate blocks of the same uorf/mrna
    geom_line(aes(y = y, group = factor(y))) +
    # draw exonic ribbon
    geom_ribbon(
      aes(ymin = ymin, ymax = ymax, group = fid, fill = type),
      show.legend = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_minimal() +  # remove unnessary plot components
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank())
}

## load CAGE dataset ################################################################################
# cage.files <- list.files('../pbio_rev/CAGE_data/tagCluster_simpleTpm',
#                          pattern = 'raw.(plus|minus).bw', full.names = TRUE)
# cage.files <- t(matrix(cage.files, nrow = 2))[, c(2, 1)]
# colnames(cage.files) <- c('cage.plus', 'cage.minus')
# rownames(cage.files) <- sub('.*/(.*?).CTSS.*', '\\1', cage.files[, 1])
# cage.files <- cage.files[c(2, 6:9,3:5, 12:15, 18:19, 10, 16, 1, 11, 17), ]

cage.files <- list.files('CAGE_bam_coverage', pattern = 'raw.(plus|minus).bw', full.names = TRUE)
cage.files <- t(matrix(cage.files, nrow = 2))[, c(2, 1)]
colnames(cage.files) <- c('cage.plus', 'cage.minus')
rownames(cage.files) <- sub('.*/(.*?).CTSS.*', '\\1', cage.files[, 1])

modencode.sra <- readxl::read_excel('./modencode_sample_name.xlsx')
setDT(modencode.sra)
modencode.sra <- modencode.sra[!is.na(name)]
modencode.sra <- modencode.sra[!duplicated(name)]
#modencode.sra[, name := make.unique(name)]

cage.files <- cage.files[modencode.sra$CAGE, ]
rownames(cage.files) <- modencode.sra$name

cage.files <- cage.files[sort(rownames(cage.files))[c(
  9, 28, 33, 1:8, 10:12, 14, 13, 15:17, 30, 29, 31:32, 18:19, 21:22, 20, 34, 24, 23, 26, 27, 25)], ]


rownames(cage.files)[rownames(cage.files) == 'pupae_WPP12hr'] <- 'white_prepupae_12hr'
rownames(cage.files)[rownames(cage.files) == 'pupae_2d_postWPP'] <- 'pupae_2days_post_white_prepupae'
rownames(cage.files)[rownames(cage.files) == 'pupae_WPP2d_fat'] <- 'pupae_2days_post_white_prepupae_fat_body'
rownames(cage.files)[rownames(cage.files) == 'pupae_WPP2d_CNS'] <- 'pupae_2days_post_white_prepupae_CNS'

## functions to make coverage plot for CAGE data ####################################################
ExtractGeneCov.cage <- function(gene, samples = rownames(cage.files)){
  # get gene range
  g.intv <- fly.genes[gene]
  chrm <- as.character(runValue(seqnames(g.intv)))
  g.strand <- as.character(runValue(strand(g.intv)))
  g.intv <- keepSeqlevels(g.intv, chrm)
  # gt bw file list
  if(g.strand == '+'){
    bw.files.list <- cage.files[samples, 1]
  }else{
    bw.files.list <- cage.files[samples, 2]
  }
  bw.files.list <- BigWigFileList(as.vector(bw.files.list))
  # extract coverage
  mat.cov <- sapply(bw.files.list, function(bw){
    bw <- import(bw, which = g.intv)
    bw.cov <- coverage(keepSeqlevels(bw, chrm), weight = 'score')
    as.vector(Views(bw.cov[[1]], ranges(g.intv))[[1]])
  })
  if(g.strand == '-'){
    mat.cov <- mat.cov[nrow(mat.cov):1, ] * (-1)
  }
  mat.cov
}

MakeCovPlot.cage <- function(g, ..., plt.samples = rownames(cage.files),
                             log.scale = FALSE, nm = FALSE, plot = TRUE){
  #g <- 'FBgn0034876'
  g.sample <- plt.samples
  g.mat <-ExtractGeneCov.cage(g, plt.samples)
  if(log.scale){
    g.mat <- apply(g.mat, 2, function(x){
      m <- min(x[x > 0])
      log(x/m + 1)
    })
  }
  colnames(g.mat) <- g.sample
  
  g.df <- as.data.frame(g.mat)
  g.df$posn <- seq_len(nrow(g.df))
  
  g.exons <- ExtractGeneInfo(g, ...)
  g.plt <- RemoveIntron(g.df, g.exons)
  
  g.plt <- melt(
    g.plt, id.vars = 'posn', variable.name = 'track', value.name = 'coverage')
  if(nm){
    setDT(g.plt)
    g.plt[, coverage := coverage / max(coverage), by = .(track)]
    setDF(g.plt)
  }
  if(is.data.table(g.plt)){
    print(head(g.plt))
  }
  #g.plt[c('type', 'sample')] <- tstrsplit(g.plt$track, '|', fixed = TRUE)
  # could also use levels = plt.samples
  g.plt$track <- factor(g.plt$track, levels = plt.samples)
  
  if(plot){
    ggplot(g.plt, aes(x = posn, y = coverage)) +
      facet_grid(track ~ ., scales = 'free_y') +
      geom_ribbon(aes(ymax = coverage, ymin = 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      theme_classic() +
      theme(strip.text.y = element_text(angle = 0),
            strip.background = element_blank(),
            axis.text = element_blank())
  }else{
    invisible(g.plt)
  }
}

MergeCageSample <- function(dt, keep.samples = NULL){
  sample.list <- list(
    em_0_2h = 'embryo_00_02hr',
    em_2_6h = c('embryo_02_04hr', 'embryo_04_06hr'),
    em_6_12h = c('embryo_06_08hr', 'embryo_08_10hr'),
    em_12_24h = c('embryo_12_14hr', 'embryo_16_18hr', 'embryo_20_22hr'),
    larva = c(
      'larva_L3_puffstage_3_6', 'larva_L3_wandering_stage_carcass', 'larva_L3_CNS',
      'larva_L3_wandering_stage_digestive_system', 'larva_L3_wandering_stage_imaginal_discs'),
    pupa = c('pupae_2days_post_white_prepupae', 'pupae_2days_post_white_prepupae_CNS',
             'pupae_2days_post_white_prepupae_fat_body'),
    female_body = c(
      'mated_female_eclosion_4_days_ovaries', "mixed_adults_eclosion_4_days_digestive_system",
      'mixed_adults_eclosion_20_days_digestive_system', 'mixed_adults_eclosion_4_days_carcass',
      'virgin_female_eclosion_4_days_ovaries'),
    male_body = c(
      'mated_male_eclosion_4_days_accessory_glands', 'mated_male_eclosion_4_days_testes',
      'mixed_adults_eclosion_20_days_digestive_system', 'mixed_adults_eclosion_4_days_carcass',
      "mixed_adults_eclosion_4_days_digestive_system"),
    female_head = c(
      'mated_female_eclosion_1_day_heads', 'mated_female_eclosion_20_days_heads'),
    male_head = c(
      'mated_male_eclosion_1_day_heads', 'mated_male_eclosion_20_days_heads'),
    S2_wd = 'S2_DRSC'
  )
  names(sample.list) <- splt.name[names(sample.list)]
  
  dt.char <- as.data.table(dt)
  dt.char[, track := as.character(track)]
  tmp <- lapply(sample.list, function(x){
    dt.char[track %in% x][, .(coverage = sum(coverage)), by = 'posn']
  })
  if(!is.null(keep.samples)){
    tmp <- tmp[keep.samples]
    res <- rbindlist(tmp, idcol = 'track')
    res[, track := factor(track, levels = keep.samples)]
  }else{
    res <- rbindlist(tmp, idcol = 'track')
    res[, track := factor(track, levels = splt.name)]
  }
  setDF(res)
  res
}
