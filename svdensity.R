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

chm13intersectfn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi.specprec.tsv"
chm13intersectvariants = read.table(chm13intersectfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
chm13variants = chm13intersectvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
goodlist <- paste('chr', c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                           '21', '22', "X"), sep = "")
binsize <- 1000000
chm13variants$STARTBIN <- floor(as.numeric(chm13variants$POS) / binsize) * binsize
groups <- chm13variants %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
unique(groups$CHROM)
groups$CHROM = factor(groups$CHROM, levels=goodlist)
groups$logn <- log2(groups$n + 1)

chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX")
chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                              12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                              21, 22, 23)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))
sizes = c(248387497,
          242696747,
          201106605,
          193575430,
          182045437,
          172126870,
          160567423,
          146259322,
          150617274,
          134758122,
          135127772,
          133324781,
          114240146,
          101219177,
          100338308,
          96330493,
          84277185,
          80542536,
          61707359,
          66210247,
          45827691,
          51353906,
          154259625)
sum(sizes)
sizedf <- data.frame(matrix(ncol = 0, nrow = 23))
sizedf$CHROM <- goodlist
sizedf$CHROM <- factor(sizedf$CHROM, levels=goodlist)
sizedf$size <- sizes
sizedf
sizedf$size

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

nonsyntenyfn <- "/home/mkirsche/JasmineT2T/chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed"
nonsynteny <- read.table(nonsyntenyfn, comment.char = "#", sep = "\t", header = F, stringsAsFactors=FALSE, colClasses = 'character')
colnames(nonsynteny) <- c("CHROM", "Start", "End")
nonsynteny <- nonsynteny %>% filter(CHROM %in% goodlist)
nonsynteny$CHROM <- factor(nonsynteny$CHROM, levels=goodlist)

plot <- ggplot() +scale_y_discrete(name = "chromosome", limits = names(chrom_key)) +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), fill = "white") +
  geom_rect(data = groups, aes(ymin = as.numeric(CHROM) - 0.3, 
                                 ymax = as.numeric(CHROM) + 0.1, 
                                 xmax = STARTBIN + binsize, xmin = STARTBIN, fill = logn)) + 
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), color = "black", fill = "NA") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "#44BB99") +
  geom_rect(data = nonsynteny, aes(xmin = as.numeric(Start), xmax = as.numeric(End), ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), fill = "#EE8866") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "NA") +
  
  scale_x_continuous(labels = scales::comma) +
  xlab("position (bp)") +
  scale_fill_gradient(low =  "white", high = "black", name = "SV Density (log2)", breaks = c(2.5, 5, 7.5, 10), limits = c(0, 10))
ggsave("/home/mkirsche/JasmineT2T/density_log.png", width = 10, height = 12)  
ggsave("/home/mkirsche/JasmineT2T/density_log.svg", width = 10, height = 12)  


ggplot() +scale_y_discrete(name = "chromosome", limits = names(chrom_key)) + theme(legend.position = "none") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), fill = "white") +
  geom_rect(data = groups, aes(ymin = as.numeric(CHROM) - 0.3, 
                               ymax = as.numeric(CHROM) + 0.1, 
                               xmax = STARTBIN + binsize, xmin = STARTBIN, fill = logn)) + 
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), color = "black", fill = "NA") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "#44BB99") +
  geom_rect(data = nonsynteny, aes(xmin = as.numeric(Start), xmax = as.numeric(End), ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), fill = "#EE8866") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "NA") +
  
  scale_x_continuous(labels = scales::comma) +
  xlab("position (bp)") +
  scale_fill_gradient(low =  "white", high = "black", name = "SV Density (log2)", limits = c(0, 10))
ggsave("/home/mkirsche/JasmineT2T/density_log_nolegend.png", width = 10, height = 12)  
ggsave("/home/mkirsche/JasmineT2T/density_log_nolegend.svg", width = 10, height = 12)  
ggsave("/home/mkirsche/JasmineT2T/density_log_nolegend.pdf", width = 10, height = 12)  

ggplot() +scale_y_discrete(name = "chromosome", limits = names(chrom_key)) +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), fill = "white") +
  geom_rect(data = groups, aes(ymin = as.numeric(CHROM) - 0.3, 
                               ymax = as.numeric(CHROM) + 0.1, 
                               xmax = STARTBIN + binsize, xmin = STARTBIN, fill = n)) + 
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), color = "black", fill = "NA") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "#44BB99") +
  geom_rect(data = nonsynteny, aes(xmin = as.numeric(Start), xmax = as.numeric(End), ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), fill = "#EE8866") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "NA") +
  
  scale_x_continuous(labels = scales::comma) +
  xlab("position (bp)") +
  scale_fill_gradient(low =  "white", high = "black", name = "SV Density")
ggsave("/home/mkirsche/JasmineT2T/density.png", width = 10, height = 12)  


# Singletons
singlegroups <- chm13variants %>% filter(SUPP=="1") %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
unique(singlegroups$CHROM)
singlegroups$CHROM = factor(singlegroups$CHROM, levels=goodlist)
singlegroups$logn <- log2(singlegroups$n + 1)

ggplot() +scale_y_discrete(name = "chromosome", limits = names(chrom_key)) + theme(legend.position = "none") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), fill = "white") +
  geom_rect(data = singlegroups, aes(ymin = as.numeric(CHROM) - 0.3, 
                               ymax = as.numeric(CHROM) + 0.1, 
                               xmax = STARTBIN + binsize, xmin = STARTBIN, fill = logn)) + 
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), color = "black", fill = "NA") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "#44BB99") +
  geom_rect(data = nonsynteny, aes(xmin = as.numeric(Start), xmax = as.numeric(End), ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), fill = "#EE8866") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "NA") +
  
  scale_x_continuous(labels = scales::comma) +
  xlab("position (bp)") +
  scale_fill_gradient(low =  "white", high = "black", name = "Singleton Density (log2)")
ggsave("/home/mkirsche/JasmineT2T/density_singletons_log_nolegend.png", width = 10, height = 12)  


ggplot() +scale_y_discrete(name = "chromosome", limits = names(chrom_key)) +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), fill = "white") +
  geom_rect(data = singlegroups, aes(ymin = as.numeric(CHROM) - 0.3, 
                               ymax = as.numeric(CHROM) + 0.1, 
                               xmax = STARTBIN + binsize, xmin = STARTBIN, fill = logn)) + 
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) - .30, ymax = as.numeric(CHROM) + .10), color = "black", fill = "NA") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "#44BB99") +
  geom_rect(data = nonsynteny, aes(xmin = as.numeric(Start), xmax = as.numeric(End), ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), fill = "#EE8866") +
  geom_rect(data = sizedf, aes(xmin = 0, xmax = size + 100000, ymin = as.numeric(CHROM) + .10, ymax = as.numeric(CHROM) + .30), color = "black", fill = "NA") +
  
  scale_x_continuous(labels = scales::comma) +
  xlab("position (bp)") +
  scale_fill_gradient(low =  "white", high = "black", name = "Singleton Density (log2)")
ggsave("/home/mkirsche/JasmineT2T/density_singletons_log.png", width = 10, height = 12)  



outfile <- "/home/mkirsche/JasmineT2T/densitylegend.svg"
#png(outfile, width = 4096, height = 4096)
Cairo(file=outfile, 
      type="svg",
      units="in", 
      width=4, 
      height=6, 
      pointsize=12, 
      dpi=72)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", 
       legend = c("Non-syntenic", "Syntenic"), 
       col = c('#EE8866', '#44BB99'), 
       pch = 15,
       bty = "n", 
       pt.cex = 3, 
       cex = 1.5, 
       text.col = "black")
dev.off()


outfile <- "/home/mkirsche/JasmineT2T/densitylegend.png"
#png(outfile, width = 4096, height = 4096)
Cairo(file=outfile, 
      type="png",
      units="in", 
      width=4, 
      height=6, 
      pointsize=12, 
      dpi=72)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", 
       legend = c("Non-syntenic", "Syntenic"), 
       col = c('#EE8866', '#44BB99'), 
       pch = 15,
       bty = "n", 
       pt.cex = 3, 
       cex = 1.5, 
       text.col = "black")
dev.off()



outfile <- "/home/mkirsche/JasmineT2T/chm13_chr_svs.png"
ggplot(groups %>% filter(CHROM %in% goodlist), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity") +
  xlab("Genomic Position") +
  ylab("SVs per 100kbp bin (log2)") +
  ggtitle("SV Density (CHM13)") +
  facet_wrap(~ CHROM, ncol = 5, scales = "free") +
  theme(plot.title = element_text(hjust = .5)) 
ggsave(outfile, width = 12, height = 12)
