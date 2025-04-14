setwd("~/syq/project/uorf/gbrowser_graph/")
source('Fig2_source.R')

sample.libsize <- readxl::read_excel('libsize.xlsx')
setDT(sample.libsize)
sample.libsize[, med := mean(size), by = .(type)]

sample.merge <- c("mature_oocyte", "em_0_2h", "em_2_6h", "em_6_12h", "em_12_24h",
                  "larva", "pupa", "female_body", "male_body", "female_head", "male_head")


## clock #########################################################################################
g <- 'FBgn0023076'
g.drop <- c("FBtr0076785", "FBtr0334648", "FBtr0100134", "FBtr0334647")
drop.manual <-uorf.annot[transcript%in%g.drop,name]

setDT(dtt.cov)
dtt.cov[, sample := tstrsplit(sample, '(?=\\d$)', perl = TRUE)[[1]] ]
dtt.cov <- dtt.cov[, .(coverage = sum(coverage)), by = .(type, sample, posn)]
dtt.cov[sample.libsize, `:=`(size = i.size, med = i.med), on = .(type, sample)]
dtt.cov[, coverage := coverage / size * med]
dtt.cov[, sample := factor(sample, levels = sample.merge)]

p1 <- ggplot(dtt.cov[type == 'mrna' & posn <= 500, ], aes(x = posn, y = coverage)) +
  facet_grid(sample ~ ., scales = 'free_y') +
  geom_ribbon(aes(ymax = coverage, ymin = 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'mRNA coverage') +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.background = element_blank(),
        axis.text = element_blank())
p2 <- GeneModelExonPlot(g, xlim = c(1, 500), drop.isoform = g.drop, drop.manual = drop.manual)
p3 <- ggplot(dtt.cov[type == 'ribo' & posn <= 500, ], aes(x = posn, y = coverage)) +
  facet_grid(sample ~ ., scales = 'free_y') +
  geom_ribbon(aes(ymax = coverage, ymin = 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'RPF coverage') +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.background = element_blank(),
        axis.text = element_blank())

m <- cowplot::plot_grid(p1, p2, p3, align = 'v', axis = 'lr', ncol = 1, rel_heights = c(4, 1, 4))
print(m)

## cycle #########################################################################################
g <- 'FBgn0023094'
g.drop <- c("")

dtt.cov <- MakeCovPlot(g, plot = FALSE, drop.isoform = g.drop,
                       plt.samples = rownames(bw.files.meta)[8:13])
setDT(dtt.cov)
dtt.cov[, sample := tstrsplit(sample, '(?=\\d$)', perl = TRUE)[[1]] ]
dtt.cov <- dtt.cov[, .(coverage = sum(coverage)), by = .(type, sample, posn)]
dtt.cov[sample.libsize, `:=`(size = i.size, med = i.med), on = .(type, sample)]
dtt.cov[, coverage := coverage / size * med]
dtt.cov[, sample := factor(sample, levels = sample.merge)]

p1 <- ggplot(dtt.cov[type == 'mrna' & posn <= 500, ], aes(x = posn, y = coverage)) +
  facet_grid(sample ~ ., scales = 'free_y') +
  geom_ribbon(aes(ymax = coverage, ymin = 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'mRNA coverage') +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.background = element_blank(),
        axis.text = element_blank())
p2 <- GeneModelExonPlot(g, xlim = c(1, 500), drop.isoform = g.drop, drop.manual = drop.manual)
p3 <- ggplot(dtt.cov[type == 'ribo' & posn <= 500, ], aes(x = posn, y = coverage)) +
  facet_grid(sample ~ ., scales = 'free_y') +
  geom_ribbon(aes(ymax = coverage, ymin = 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = 'RPF coverage') +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0, size = 10),
        strip.background = element_blank(),
        axis.text = element_blank())

m <- cowplot::plot_grid(p1, p2, p3, align = 'v', axis = 'lr', ncol = 1, rel_heights = c(4, 1, 4))
print(m)
