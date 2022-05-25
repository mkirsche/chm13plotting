library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(reshape2)
library(gridExtra)
library(rlist)
library(VennDiagram)
library(Cairo)
library(ggpubr)
library(here)
library(svglite)

make_pca_scatter <- function(variants, pref, samples)
{
  variants$SEP <- sub("\\s+$", "", gsub('(.{1})', '\\1 ', variants$SUPP_VEC))
  samples$SAMPLE
  sepdf <- variants %>% separate(SEP, unlist(samples$SAMPLE), sep=" ")
  pcadf <- sepdf[,unlist(samples$SAMPLE)]
  colnames(pcadf)
  pcamatrix<-t(data.matrix(pcadf))
  pca <- prcomp(pcamatrix, center = TRUE)
  df_out <- as.data.frame(pca$x)
  df_out$SAMPLE <- rownames(df_out)
  df_out <- merge(df_out, samples, by = "SAMPLE")
  ggplot(df_out, aes(x = PC1, y = PC2, label = SAMPLE, color = POPULATION)) + geom_point() + geom_text(nudge_y = -1.5, size=2)
  outfile <- paste("/home/mkirsche/JasmineT2T/", pref, "_", "pcapop.png", sep = "")
  ggsave(outfile, width = 8, height = 8)
  ggplot(df_out, aes(x = PC1, y = PC2, label = SAMPLE, color = SUPERPOPULATION)) + geom_point() + geom_text(nudge_y = -1.5, size=2)
  outfile <- paste("/home/mkirsche/JasmineT2T/", pref, "_", "pcasuperpop.png", sep = "")
  ggsave(outfile, width = 8, height = 8)
  
  # Get the number of SV calls in each sample
  sampletally <- sepdf %>% gather(x, value, unlist(samples$SAMPLE)) %>% group_by(x) %>% filter(value==1)
  names(sampletally)[names(sampletally) == 'x'] <- 'SAMPLE'
  nrow(sampletally)
  colnames(sampletally)
  colnames(samples)
  summarized <- sampletally %>% group_by(SVTYPE, SAMPLE) %>% summarise(counts=n())
  nrow(summarized)
  summarized
  summarized <- merge(summarized, samples, by = "SAMPLE")
  summarized
  sampletally$TECH
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  ggplot(summarized, aes(x = paste(SAMPLE, paste("(",POPULATION, ")",sep=""), sep="\n"), fill = SVTYPE, y=counts)) + geom_bar(stat = "identity") + xlab("Sample") + ylab("Number of Variants") + facet_grid(cols=vars(SUPERPOPULATION), scales = "free_x", space="free_x") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 8, angle = 30),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          #legend.text = element_blank(),
          text = element_text(size = 16),
          #legend.position = "none"
    ) + 
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    guides(fill=guide_legend(title="Type"))
  outfile <- paste("/home/mkirsche/JasmineT2T/", pref, "_", "samplecountsbar.png", sep = "")
  outfile <- paste("/home/mkirsche/JasmineT2T/", pref, "_", "samplecountsbar.svg", sep = "")
  
  ggsave(outfile, width = 12, height = 8)
  
  
  
  sampletally <- sepdf %>% filter(abs(as.numeric(SVLEN)) >= 50 | SVTYPE == "TRA") %>% gather(x, value, unlist(samples$SAMPLE)) %>% group_by(x) %>% filter(value==1)
  names(sampletally)[names(sampletally) == 'x'] <- 'SAMPLE'
  nrow(sampletally)
  colnames(sampletally)
  colnames(samples)
  summarized <- sampletally %>% group_by(SVTYPE, SAMPLE) %>% summarise(counts=n())
  nrow(summarized)
  summarized
  summarized <- merge(summarized, samples, by = "SAMPLE")
  summarized
  sampletally$TECH
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  ggplot(summarized, aes(x = paste(SAMPLE, paste("(",POPULATION, ")",sep=""), sep="\n"), fill = SVTYPE, y=counts)) + geom_bar(stat = "identity") + xlab("Sample") + ylab("Number of Variants") + facet_grid(cols=vars(SUPERPOPULATION), scales = "free_x", space="free_x") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 8, angle = 30),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          #legend.text = element_blank(),
          text = element_text(size = 16),
          #legend.position = "none"
    ) + 
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    guides(fill=guide_legend(title="Type"))
  outfile <- paste("/home/mkirsche/JasmineT2T/", pref, "_", "longsamplecountsbar.png", sep = "")
  outfile <- paste("/home/mkirsche/JasmineT2T/", pref, "_", "longsamplecountsbar.svg", sep = "")
  
  ggsave(outfile, width = 12, height = 8)
}

fn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi.specprec.tsv"
samplesfn <- "/home/mkirsche/JasmineT2T/filelist.txt"

variants <-  read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
samples <- read.table(samplesfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
pref <- "chm13_intersect_hifi"
make_pca_scatter(variants, pref, samples)
