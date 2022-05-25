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
library(ggpattern)

chm13intersectannofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_anno_annotated.tsv"
df = read.table(chm13intersectannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#plot_length_line_highlight(chm13intersectannovariants, "CHM13", highlighted_length_line_ofn, T, .92, .85)

df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
summarized <- df %>% group_by(SVTYPE, INTERSECTS_NONSYNTENY_CHM13) %>% summarise(counts=n()) %>% filter(SVTYPE == "INS" | SVTYPE == "DEL")
summarized$CATEGORY <- paste(summarized$SVTYPE, "_", summarized$INTERSECTS_NONSYNTENY_CHM13, sep = "")
summarized$CATEGORY <- ifelse(summarized$CATEGORY == "INS_0", "Syntenic",summarized$CATEGORY)
summarized$CATEGORY <- ifelse(summarized$CATEGORY == "INS_1", "Non-syntenic",summarized$CATEGORY)
summarized$CATEGORY <- ifelse(summarized$CATEGORY == "DEL_0", "Syntenic",summarized$CATEGORY)
summarized$CATEGORY <- ifelse(summarized$CATEGORY == "DEL_1", "Non-syntenic",summarized$CATEGORY)
summarized
ggplot(data = summarized, aes(x = SVTYPE, y = counts, fill = SVTYPE, pattern = CATEGORY)) +
  geom_bar_pattern(stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  #scale_fill_manual(values = c('#DF4828', '#4EB265')) +
  scale_fill_manual(values = c('#CC6677', '#44AA99')) +
  scale_pattern_manual(values = c("stripe", "none")) +
  theme_classic() +
  ylab("Count") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        #legend.text = element_blank(),
        legend.title = element_blank(),
        #legend.position = "none",
  ) +
  #ylim(0, 70000) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern_newcolors.pdf"
ggsave(ofn, width= 4, height = 8)
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern_newcolors.png"
ggsave(ofn, width= 4, height = 8)
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern_newcolors.svg"
ggsave(ofn, width= 4, height = 8)

ggplot(data = summarized, aes(x = SVTYPE, y = counts, fill = SVTYPE, pattern = CATEGORY)) +
  geom_bar_pattern(stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = c('#DF4828', '#4EB265')) +
  scale_pattern_manual(values = c("stripe", "none")) +
  theme_classic() +
  ylab("Count") +
  theme(plot.title = element_text(size = 24, hjust = 0.5),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        #legend.text = element_blank(),
        #legend.title = element_blank(),
        legend.position = "none",
  ) +
  #ylim(0, 70000) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern_nolegend.pdf"
ggsave(ofn, width= 3, height = 8)
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern_nolegend.png"
ggsave(ofn, width= 3, height = 8)
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern_nolegend.svg"
ggsave(ofn, width= 3, height = 8)



plot_table <- function(df, title)
{
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  df$LenCategory = ifelse(df$ABSLEN >= 0 & df$ABSLEN <= 50, "30-50", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 50 & df$ABSLEN <= 100, "50-100", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 100 & df$ABSLEN <= 200, "100-200", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 200 & df$ABSLEN <= 300, "200-300", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 300 & df$ABSLEN <= 400, "300-400", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 400 & df$ABSLEN <= 500, "400-500", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 500 & df$ABSLEN <= 750, "500-750", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 750 & df$ABSLEN <= 1000, "750-1k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 1000 & df$ABSLEN <= 2000, "1k-2k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 2000 & df$ABSLEN <= 5000, "2k-5k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 5000 & df$ABSLEN <= 10000, "5k-10k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 10000, "10k+", df$LenCategory)
  
  df$LenCategory = ifelse(df$SVTYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  plot <- ggplot(summarized) +
    geom_bar(data = summarized %>% filter(SVTYPE != "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = counts))+
    geom_bar(data = summarized %>% filter(SVTYPE == "DEL"), position = "stack", stat = "identity", aes(x = LenCategory, fill = SVTYPE, y = -counts))+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    labs(title = paste("Indels by Length (", title, ")", sep = "")) +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 14, angle=30),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          text = element_text(size = 16),
          legend.position = "None",
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  return(plot)
}

plot_allele_frequencies <- function(variants, outfile, title, maxy) {
  variants$SUPP_INT = strtoi(variants$SUPP)
  
  childcounts <- variants %>% filter(variants$SUPP_INT == 17) %>% group_by(SUPP_INT) %>% count()
  childcounts$SVTYPE = "INS"
  childcounts$shape = 23
  
  ggplot(variants, aes(x = SUPP_INT, fill = SUPP_INT)) +
    geom_bar(position = "stack", stat = "count", fill = "#2093c3") +
    ggtitle(title) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 24)) +
    xlab("Allele Frequency") +
    ylab("Number of Variants")
    variants$SVTYPE <- ifelse(variants$SVTYPE == "INS" | variants$SVTYPE == "DEL", variants$SVTYPE, "Other")
    ggsave(outfile, height = 8, width = 8)
    colorpalette <- c(DEL = '#CC6677',Other = '#1965B0',INS = '#44AA99')
  #colorpalette <- c(DEL = '#DF4828',Other = '#1965B0',INS = '#4EB265')
  variants$SUPP_INT = strtoi(variants$SUPP)
  ggplot(variants, aes(x = SUPP_INT, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "count",) +
    ggtitle(title) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.position = c(.9, .89),
          plot.title = element_text(hjust = 0.5, size = 24)) +
    xlab("Allele Frequency") +
    ylab("Number of Variants") +
    geom_point(data = childcounts, aes(x = SUPP_INT, y=n + 1500), shape = 25, fill = "darkred", color = "darkred", size = 5) +
    ylim(0, maxy) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette)
    
    ggsave(paste(outfile, "_type.png", sep = ""), height = 8, width = 10)
    ggsave(paste(outfile, "_type.svg", sep = ""), height = 8, width = 10)
}

suppvec_hist_highlightdiscordant2 <- function(df, caller,div) {
  offset <- nrow(df)/div
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", "Child Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", "Father Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", "Mother Only", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", "Son/Father", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", "Son/Mother", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", "Both parents", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- factor(df$SUPP_VEC_STRING,levels = c("All three", "Both parents", "Father Only", 
                                                             "Mother Only", "Son/Father", "Son/Mother", "Child Only"))
  #colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  colorpalette <- c(DEL = '#DF4828', DUP = '#1965B0',INS = '#4EB265', TRA = '#7BAFDE', INV = '#F7CB45')
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  
  childcounts <- df %>% filter(df$SUPP_VEC_STRING == "Child Only") %>% group_by(SUPP_VEC_STRING) %>% count()
  childcounts$SVTYPE = "INS"
  childcounts$shape = 23
  
  summarized <- df %>% group_by(SVTYPE, SUPP_VEC_STRING) %>% summarise(counts=n())
  
  
  plot <- ggplot(summarized, aes(x = SUPP_VEC_STRING, y = counts, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "identity") +
    labs(title = paste("SVs by Sample Presence (", caller, ")", sep = "")) +
    xlab("Samples") +
    ylab("Count") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(plot.title = element_text(size = 18, hjust = 0.5),
          axis.text.x = element_text(size = 12, angle = 30, margin = margin(t = 15)),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 16),
          legend.position = "None"#c(legendx, legendy),
    ) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_point(data = childcounts, aes(x = SUPP_VEC_STRING, y=n+offset), shape = 25, fill = "darkred", color = "darkred", size = 5)
  return(plot)
}

suppvec_hist_highlightdiscordant <- function(df, caller, outfile, filter, offset, legendx, legendy, title, titlesize, ylim, child, parent1, parent2) {
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df$SUPP_VEC_STRING <- str_pad(as.character(df$SUPP_VEC), 3, "left", "0")
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "100", paste(child, " Only", sep = ""), df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "010", paste(parent1, " Only", sep = ""), df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "001", paste(parent2, " Only", sep = ""), df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "110", paste(child, "/", parent1, sep = ""), df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "101", paste(child, "/", parent2, sep = ""), df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "011", paste(parent1, "/", parent2, sep = ""), df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- ifelse(df$SUPP_VEC_STRING == "111", "All three", df$SUPP_VEC_STRING)
  df$SUPP_VEC_STRING <- factor(df$SUPP_VEC_STRING,levels = c("All three", paste(parent1, "/", parent2, sep = ""), paste(parent1, " Only", sep = ""), 
                                                             paste(parent2, " Only", sep = ""), paste(child, "/", parent1, sep = ""),  paste(child, "/", parent2, sep = ""), paste(child, " Only", sep = "")))
  
  suppveccounts <- df %>% count(SUPP_VEC_STRING)
  suppveccounts
  suppveccounts$SVTYPE = "INS"
  colorpalette <- c(DEL = '#DF4828', DUP = '#1965B0',INS = '#4EB265', TRA = '#7BAFDE', INV = '#F7CB45')
  #colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  childcounts <- df %>% filter(df$SUPP_VEC_STRING == paste(child, " Only", sep = "")) %>% group_by(SUPP_VEC_STRING) %>% count()
  childcounts$SVTYPE = "INS"
  childcounts$shape = 23
  ggplot(df, aes(x = SUPP_VEC_STRING, fill = SVTYPE)) +
    geom_bar(position = "stack", stat = "count") +
    xlab("Samples") +
    ylab("Count") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    theme(axis.text.x = element_text(size = 16, angle = 30, margin = margin(t = 21)),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          legend.position = c(legendx, legendy),
          plot.title = element_text(size = titlesize, hjust = 0.5),
    ) +
    ggtitle(paste("Mendelian Discordance (", title, ")", sep = "")) +
    ylim(0, ylim) +
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = suppveccounts, aes(x = SUPP_VEC_STRING, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_point(data = childcounts, aes(x = SUPP_VEC_STRING, y=n+offset), shape = 25, fill = "darkred", color = "darkred", size = 5)
  ggsave(outfile, width= 9, height = 8)
}

plot_length <- function(df, caller, outfile, filter, legendx, legendy) 
{
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  df$LenCategory = ifelse(df$ABSLEN >= 0 & df$ABSLEN <= 50, "30-50", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 50 & df$ABSLEN <= 100, "50-100", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 100 & df$ABSLEN <= 200, "100-200", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 200 & df$ABSLEN <= 300, "200-300", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 300 & df$ABSLEN <= 400, "300-400", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 400 & df$ABSLEN <= 500, "400-500", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 500 & df$ABSLEN <= 750, "500-750", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 750 & df$ABSLEN <= 1000, "750-1k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 1000 & df$ABSLEN <= 2000, "1k-2k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 2000 & df$ABSLEN <= 5000, "2k-5k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 5000 & df$ABSLEN <= 10000, "5k-10k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 10000, "10k+", df$LenCategory)
  
  df$LenCategory = ifelse(df$SVTYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  colorpalette <- c(DEL = '#DF4828', DUP = '#1965B0',INS = '#4EB265', TRA = '#7BAFDE', INV = '#F7CB45')
  #colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  
  ggplot(summarized, aes(fill = SVTYPE, x = LenCategory)) +
    geom_bar(data = summarized %>% filter(SVTYPE != "DEL"), stat = "identity", aes(y = counts), position = "stack")+
    geom_bar(data = summarized %>% filter(SVTYPE == "DEL"), stat = "identity", aes(y = -counts), position = "stack")+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    ggtitle(paste('SV Size Distribution (', caller, ')', sep = '')) +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_fill_manual(name = "SVTYPE", values = colorpalette) +
    geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
}

#df<-chm13intersectvariants
#caller <- "CHM13"
#outfile <- chm13_intersect_length_line_ofn
#filter <- T
#legendx <- .92
#legendy <- .85
plot_length_line <- function(df, caller, outfile, filter, legendx, legendy) 
{
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  df$LenCategory = ifelse(df$ABSLEN >= 0 & df$ABSLEN <= 50, "30-50", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 50 & df$ABSLEN <= 100, "50-100", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 100 & df$ABSLEN <= 200, "100-200", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 200 & df$ABSLEN <= 300, "200-300", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 300 & df$ABSLEN <= 400, "300-400", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 400 & df$ABSLEN <= 500, "400-500", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 500 & df$ABSLEN <= 750, "500-750", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 750 & df$ABSLEN <= 1000, "750-1k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 1000 & df$ABSLEN <= 2000, "1k-2k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 2000 & df$ABSLEN <= 5000, "2k-5k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 5000 & df$ABSLEN <= 10000, "5k-10k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 10000, "10k+", df$LenCategory)
  
  df$LenCategory = ifelse(df$SVTYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  
  colorpalette <- c(DEL = '#DF4828', DUP = '#1965B0',INS = '#4EB265', TRA = '#7BAFDE', INV = '#F7CB45')
  #colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  summarized$LenCategory <- factor(summarized$LenCategory,levels = labellist)
  
  summarized
  lineplot <- ggplot(summarized %>% filter(SVTYPE == "INS" | SVTYPE == "DEL"), aes(fill = SVTYPE, x = LenCategory, group = SVTYPE, y = counts)) +
    geom_line(aes(color = SVTYPE), size=2)+
    geom_point(aes(color = SVTYPE), size=4)+
    #geom_line(data = summarized %>% filter(SVTYPE == "DEL"), aes(y = counts), group = 2, color = '#CC6677')+
    #geom_point(data = summarized %>% filter(SVTYPE == "DEL"), color = '#CC6677')+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    ylim(0, 30000) +
    ggtitle(paste('SV Indel Balance (', caller, ')', sep = '')) +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_color_manual(values = c('#CC6677', '#44AA99')) #+
    #scale_color_manual(values = c('#DF4828', '#4EB265')) #+
    #geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
    #geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 8, height = 8)
  summarized <- df %>% group_by(SVTYPE) %>% summarise(counts=n()) %>% filter(SVTYPE == "INS" | SVTYPE == "DEL")
  summarized
  barplot <- 
    ggplot(data = summarized, aes(x = SVTYPE, y = as.numeric(counts), fill = SVTYPE)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c('#CC6677', '#44AA99')) +
    #scale_fill_manual(values = c('#DF4828', '#4EB265')) +
    theme_classic() +
    ylab("Count") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          legend.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
    ) +
    ylim(0, 70000)
  ggsave(paste(outfile, "_bar.png", sep = ""), width= 3, height = 8)
  ggsave(paste(outfile, "_bar.svg", sep = ""), width= 3, height = 8)
  
  ggarrange(barplot, lineplot, nrow = 1, ncol = 2, widths = c (1, 4))
  ggsave(paste(outfile, "_full.png", sep = ""), width= 12, height = 8)
  ggsave(paste(outfile, "_full.svg", sep = ""), width= 12, height = 8)
}

#df <- chm13annovariants
#caller <- "CHM13"
#outfile <- highlighted_length_line_ofn
#filter <- T
plot_length_line_highlight <- function(df, caller, outfile, filter, legendx, legendy) 
{
  if (filter)
  {
    df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  df$LenCategory = ifelse(df$ABSLEN >= 0 & df$ABSLEN <= 50, "30-50", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 50 & df$ABSLEN <= 100, "50-100", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 100 & df$ABSLEN <= 200, "100-200", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 200 & df$ABSLEN <= 300, "200-300", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 300 & df$ABSLEN <= 400, "300-400", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 400 & df$ABSLEN <= 500, "400-500", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 500 & df$ABSLEN <= 750, "500-750", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 750 & df$ABSLEN <= 1000, "750-1k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 1000 & df$ABSLEN <= 2000, "1k-2k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 2000 & df$ABSLEN <= 5000, "2k-5k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 5000 & df$ABSLEN <= 10000, "5k-10k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 10000, "10k+", df$LenCategory)
  
  df$LenCategory = ifelse(df$SVTYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory) %>% count()
  delcounts$SVTYPE = "DEL"
  
  nondelcounts <- df %>% filter(df$SVTYPE != "DEL") %>% group_by(LenCategory) %>% count()
  nondelcounts$SVTYPE = "INS"
  
  colorpalette <- c(DEL = '#DF4828', DUP = '#1965B0',INS = '#4EB265', TRA = '#7BAFDE', INV = '#F7CB45')
  #colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  summarized <- df %>% group_by(SVTYPE, LenCategory, INTERSECTS_NONSYNTENY_CHM13) %>% summarise(counts=n())
  summarized$LenCategory <- factor(summarized$LenCategory,levels = labellist)
  
  summarized
  lineplot <- ggplot(summarized %>% filter(SVTYPE == "INS" | SVTYPE == "DEL"), aes(fill = SVTYPE, x = LenCategory, group = SVTYPE, y = counts)) +
    geom_line(aes(color = SVTYPE), size=2)+
    geom_point(aes(color = SVTYPE), size=4)+
    #geom_line(data = summarized %>% filter(SVTYPE == "DEL"), aes(y = counts), group = 2, color = '#CC6677')+
    #geom_point(data = summarized %>% filter(SVTYPE == "DEL"), color = '#CC6677')+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("Count") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    ylim(0, 30000) +
    ggtitle(paste('SV Indel Balance (', caller, ')', sep = '')) +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_color_manual(values = c('#DF4828', '#4EB265')) #+
  #geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
  #geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
  summarized <- df %>% group_by(SVTYPE, INTERSECTS_NONSYNTENY_CHM13) %>% summarise(counts=n()) %>% filter(SVTYPE == "INS" | SVTYPE == "DEL")
  summarized$CATEGORY <- paste(summarized$SVTYPE, "_", summarized$INTERSECTS_NONSYNTENY_CHM13, sep = "")
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "INS_0", "Syntenic Insertion",summarized$CATEGORY)
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "INS_1", "Non-syntenic Insertion",summarized$CATEGORY)
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "DEL_0", "Syntenic Deletion",summarized$CATEGORY)
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "DEL_1", "Non-syntenic Deletion",summarized$CATEGORY)
  summarized
  barplot <- 
    ggplot(data = summarized %>% filter(SVTYPE == "INS" | SVTYPE == "DEL"), aes(x = SVTYPE, y = counts, fill = CATEGORY)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c('#FFCCCC', '#CCDDAA', '#CC6677', '#44AA99')) +
    theme_classic() +
    ylab("Count") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          #legend.text = element_blank(),
          #legend.title = element_blank(),
          legend.position = "none",
    ) +
    ylim(0, 70000)
  ggsave(paste(outfile, "_bar.png", sep = ""), width= 3, height = 8)
  ggsave(paste(outfile, "_bar.svg", sep = ""), width= 3, height = 8)
  
  ggarrange(barplot, lineplot, nrow = 1, ncol = 2, widths = c (1, 4))
  ggsave(paste(outfile, "_full.png", sep = ""), width= 12, height = 8)
  ggsave(paste(outfile, "_full.svg", sep = ""), width= 12, height = 8)
  
}

#df1 <- grch38intersectvariantsfiltered
#df2 <- chm13intersectvariantsfiltered
#caller <- "wm+mm2"
#outfile <- "/home/mkirsche/JasmineT2T/test.png"
compare_indel_balance <- function(df1, df2, caller, outfile, filter, legendx, legendy) 
{
  if (filter)
  {
    df1 <- df1 %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
    df2 <- df2 %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  }
  df1$ref <- "GRCh38"
  df2$ref <- "CHM13"
  df <- rbind(df1 %>% select(SVLEN, SVTYPE, ref), df2 %>% select(SVLEN, SVTYPE, ref))
  df$LenCategory = "TRA"
  
  df$ABSLEN = abs(as.numeric(df$SVLEN))
  
  df$LenCategory = ifelse(df$ABSLEN >= 0 & df$ABSLEN <= 50, "30-50", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 50 & df$ABSLEN <= 100, "50-100", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 100 & df$ABSLEN <= 200, "100-200", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 200 & df$ABSLEN <= 300, "200-300", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 300 & df$ABSLEN <= 400, "300-400", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 400 & df$ABSLEN <= 500, "400-500", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 500 & df$ABSLEN <= 750, "500-750", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 750 & df$ABSLEN <= 1000, "750-1k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 1000 & df$ABSLEN <= 2000, "1k-2k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 2000 & df$ABSLEN <= 5000, "2k-5k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 5000 & df$ABSLEN <= 10000, "5k-10k", df$LenCategory)
  df$LenCategory = ifelse(df$ABSLEN > 10000, "10k+", df$LenCategory)
  
  df$LenCategory = ifelse(df$SVTYPE == "TRA", "TRA", df$LenCategory)
  
  labellist = c("30-50", "50-100", "100-200", "200-300", "300-400", "400-500", "500-750", "750-1k", "1k-2k", "2k-5k", "5k-10k", "10k+", "TRA")
  df$LenCategory <- factor(df$LenCategory,levels = labellist)
  
  df %>% group_by(df$LenCategory) %>% tally()
  
  delcounts <- df %>% filter(df$SVTYPE == "DEL") %>% group_by(LenCategory, ref) %>% count()
  delcounts$SVTYPE = "DEL"
  
  inscounts <- df %>% filter(df$SVTYPE == "INS") %>% group_by(LenCategory, ref) %>% count()
  inscounts$SVTYPE = "INS"
  
  overall <- merge(inscounts, delcounts, by = c("LenCategory", "ref"))
  overall
  
  colorpalette <- c(DEL = '#DF4828', DUP = '#1965B0',INS = '#4EB265', TRA = '#7BAFDE', INV = '#F7CB45')
  #colorpalette <- c(DEL = '#CC6677',DUP = '#332288',INS = '#44AA99',TRA = '#88CCEE',INV = '#DDCC77')
  
  #summarized <- df %>% group_by(SVTYPE, LenCategory) %>% summarise(counts=n())
  #summarized$LenCategory <- factor(summarized$LenCategory,levels = labellist)
  
  #summarized
  ggplot(overall, aes(fill = ref, x = LenCategory, group = ref, y = (n.x - n.y))) +
    geom_line(aes(color = ref))+
    geom_point(aes(color = ref))+
    #geom_line(data = summarized %>% filter(SVTYPE == "DEL"), aes(y = counts), group = 2, color = '#CC6677')+
    #geom_point(data = summarized %>% filter(SVTYPE == "DEL"), color = '#CC6677')+
    scale_x_discrete(labels=labellist) +
    xlab("Length") +
    ylab("#INS - #DEL") +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.text.x = element_text(size = 18, angle = 30, margin = margin(t = 17, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          text = element_text(size = 20),
          legend.position = c(legendx, legendy),
    ) +
    ggtitle(paste('SV Indel Balance (', caller, ')', sep = '')) +
    #scale_fill_brewer(name = "Type", palette = "Set2") + 
    scale_color_manual(values = c('#332288', 'darkorange')) #+
  #geom_text(data = nondelcounts, aes(x = LenCategory, y=n, label=n), position=position_dodge(width=0.9), vjust=-0.75) +
  #geom_text(data = delcounts, aes(x = LenCategory, y=-1*n, label=n), position=position_dodge(width=0.9), vjust=1.5)
  
  ggsave(outfile, width= 10, height = 8)
}

###
### INPUT - GRCH38/WINNOWMAP/HIFI
###

#grch38fn <- "/home/mkirsche/JasmineT2T/grch38_wm_hifi.specprec.tsv"
#grch38variants = read.table(grch38fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38variantsfiltered <- grch38variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(grch38variantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(grch38variantsfiltered %>% filter(SVTYPE == "DEL"))

###
### INPUT - GRCH38/NEWWINNOWMAP/HIFI
###

#grch38newfn <- "/home/mkirsche/JasmineT2T/newwm/grch38_wm_hifi.specprec.tsv"
#grch38newvariants = read.table(grch38newfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38newvariantsfiltered <- grch38newvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(grch38newvariantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(grch38newvariantsfiltered %>% filter(SVTYPE == "DEL"))

###
### INPUT - GRCH38/MIX/HIFI
###

grch38intersectfn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi.specprec.tsv"
grch38intersectvariants = read.table(grch38intersectfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
grch38intersectvariantsfiltered <- grch38intersectvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(grch38intersectvariantsfiltered)
nrow(grch38intersectvariantsfiltered %>% filter(SVTYPE == "DEL"))
nrow(grch38intersectvariantsfiltered %>% filter(SVTYPE == "INS" & abs(as.numeric(SVLEN)) >= 50))
nrow(grch38intersectvariantsfiltered %>% filter(SVTYPE == "DEL" & abs(as.numeric(SVLEN)) >= 50))


###
### INPUT - CHM13/MIX/HIFI
###

chm13intersectfn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi.specprec.tsv"
chm13intersectvariants = read.table(chm13intersectfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
chm13intersectvariantsfiltered <- chm13intersectvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(chm13intersectvariantsfiltered)
nrow(chm13intersectvariantsfiltered %>% filter(SVTYPE == "DEL"))
nrow(chm13intersectvariantsfiltered %>% filter(SVTYPE == "INS" & abs(as.numeric(SVLEN)) >= 50))
nrow(chm13intersectvariantsfiltered %>% filter(SVTYPE == "DEL" & abs(as.numeric(SVLEN)) >= 50))

nrow(grch38intersectvariantsfiltered %>% filter(SUPP == "17"))
nrow(chm13intersectvariantsfiltered %>% filter(SUPP == "17"))


###
### INPUT - GRCH38/ASH/CROSSTECH
###

#grch38ashctfn <- "/home/mkirsche/JasmineT2T/grch38_ash_crosstech.tsv"
#grch38ashctvariants = read.table(grch38ashctfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38ashctvariantsfiltered <- grch38ashctvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(grch38ashctvariantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(grch38ashctvariantsfiltered %>% filter(SVTYPE == "DEL"))
#grch38_ash_ct_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_ash_ct_suppvec.png"
#suppvec_hist_highlightdiscordant(grch38ashctvariantsfiltered, 'GRCh38 Crosstech', grch38_ash_ct_suppvec_ofn, T, 2000, .9, .85, 'GRCh38 HG002', 20)


###
### INPUT - CHM13/ASH/CROSSTECH
###

#chm13ashctfn <- "/home/mkirsche/JasmineT2T/chm13_ash_crosstech.tsv"
#chm13ashctvariants = read.table(chm13ashctfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13ashctvariantsfiltered <- chm13ashctvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(chm13ashctvariantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(chm13ashctvariantsfiltered %>% filter(SVTYPE == "DEL"))
#chm13_ash_ct_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ash_ct_suppvec.png"
#suppvec_hist_highlightdiscordant(chm13ashctvariantsfiltered, 'CHM13 Crosstech', chm13_ash_ct_suppvec_ofn, T, 2000, .9, .85, 'CHM13 HG002', 20)


###
### INPUT - GRCH38/HAN/CROSSTECH
###

#grch38hanctfn <- "/home/mkirsche/JasmineT2T/grch38_han_crosstech.tsv"
#grch38hanctvariants = read.table(grch38hanctfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38hanctvariantsfiltered <- grch38hanctvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(grch38hanctvariantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(grch38hanctvariantsfiltered %>% filter(SVTYPE == "DEL"))
#grch38_han_ct_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_han_ct_suppvec.png"
#suppvec_hist_highlightdiscordant(grch38hanctvariantsfiltered, 'GRCh38 Crosstech', grch38_han_ct_suppvec_ofn, T, 2000, .9, .85, 'GRCh38 HG005', 20)


###
### INPUT - CHM13/HAN/CROSSTECH
###

#chm13hanctfn <- "/home/mkirsche/JasmineT2T/chm13_han_crosstech.tsv"
#chm13hanctvariants = read.table(chm13hanctfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13hanctvariantsfiltered <- chm13hanctvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(chm13hanctvariantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(chm13hanctvariantsfiltered %>% filter(SVTYPE == "DEL"))
#chm13_han_ct_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_han_ct_suppvec.png"
#suppvec_hist_highlightdiscordant(chm13hanctvariantsfiltered, 'CHM13 Crosstech', chm13_han_ct_suppvec_ofn, T, 2000, .9, .85, 'CHM13 HG005', 20)


###
### INDEL BALANCE - GRCH38/WINNOWMAP/HIFI
###

#grch38_length_ofn <- "/home/mkirsche/JasmineT2T/grch38_wm_hifi.lengths.png"
#plot_length(grch38variants, "GRCh38", grch38_length_ofn, T, .92, .85)
#grch38_length_line_ofn <- "/home/mkirsche/JasmineT2T/grch38_wm_hifi.lengths_line.png"
#plot_length_line(grch38variants, "GRCh38", grch38_length_line_ofn, T, .92, .85)


###
### ALLELE FREQUENCIES - GRCH38/WINNOWMAP/HIFI
###

#grch38_af_ofn <- "/home/mkirsche/JasmineT2T/grch38_wm_hifi.allelefrequency.png"
#plot_allele_frequencies(grch38variants, grch38_af_ofn, "SV Allele Frequencies (GRCh38, HiFi, wm)")


###
### INDEL BALANCE - GRCH38/MIX/HIFI
###

#grch38_intersect_length_ofn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi.lengths.png"
#plot_length(grch38intersectvariants, "GRCh38, wm+mm2", grch38_intersect_length_ofn, T, .92, .85)
grch38_intersect_length_line_ofn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi.lengths_line.png"
plot_length_line(grch38intersectvariants, "GRCh38", grch38_intersect_length_line_ofn, T, .92, .85)
grch38_intersect_length_line_ofn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi.lengths_line.pdf"
plot_length_line(grch38intersectvariants, "GRCh38", grch38_intersect_length_line_ofn, T, .92, .85)

grch38_intersect_length_line_ofn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi.lengths_line_newcolors.png"
plot_length_line(grch38intersectvariants, "GRCh38", grch38_intersect_length_line_ofn, T, .92, .85)

###
### ALLELE FREQUENCIES - GRCH38/MIX/HIFI
###

grch38_intersect_af_ofn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi_newcolors.allelefrequency.png"
plot_allele_frequencies(grch38intersectvariants, grch38_intersect_af_ofn, "SV Allele Frequencies (GRCh38)", 40000)


###
### INDEL BALANCE - CHM13/MIX/HIFI
###

#chm13_intersect_length_ofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi.lengths.png"
#plot_length(chm13intersectvariants, "CHM13, wm+mm2", chm13_intersect_length_ofn, T, .92, .85)
chm13_intersect_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi.lengths_line.png"
plot_length_line(chm13intersectvariants, "CHM13", chm13_intersect_length_line_ofn, T, .92, .85)
chm13_intersect_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi.lengths_line.pdf"
plot_length_line(chm13intersectvariants, "CHM13", chm13_intersect_length_line_ofn, T, .92, .85)
chm13_intersect_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi.lengths_line.pdf"
plot_length_line(chm13intersectvariants, "CHM13", chm13_intersect_length_line_ofn, T, .92, .85)

###
### ALLELE FREQUENCIES - CHM13/MIX/HIFI
###

chm13_intersect_af_ofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi_newcolors.allelefrequency.png"
plot_allele_frequencies(chm13intersectvariants, chm13_intersect_af_ofn, "SV Allele Frequencies (CHM13)", 40000)


###
### INPUT - GRCH38/MINIMAP2/HIFI
###

#grch38mm2fn <- "/home/mkirsche/JasmineT2T/grch38_mm2_hifi.merged.tsv"
#grch38mm2variants = read.table(grch38mm2fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38mm2variantsfiltered <- grch38mm2variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(grch38mm2variantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(grch38mm2variantsfiltered %>% filter(SVTYPE == "DEL"))


###
### INDEL BALANCE - GRCH38/MINIMAP2/HIFI
###

#grch38_mm2_length_ofn <- "/home/mkirsche/JasmineT2T/grch38_mm2_hifi.lengths.png"
#plot_length(grch38mm2variants, "GRCh38, mm2", grch38_mm2_length_ofn, T, .92, .85)
#grch38_mm2_length_line_ofn <- "/home/mkirsche/JasmineT2T/grch38_mm2_hifi.lengths_line.png"
#plot_length_line(grch38mm2variants, "GRCh38, mm2", grch38_mm2_length_line_ofn, T, .92, .85)


###
### ALLELE FREQUENCIES - GRCH38/MINIMAP2/HIFI
###

#grch38_mm2_af_ofn <- "/home/mkirsche/JasmineT2T/grch38_mm2_hifi.allelefrequency.png"
#plot_allele_frequencies(grch38mm2variants, grch38_mm2_af_ofn, "SV Allele Frequencies (GRCh38, HiFi, mm2)")


###
### INDEL BALANCE HG002 - GRCH38/WINNOWMAP/HIFI
###

#grch38hg002fn <- "/home/mkirsche/JasmineT2T/HG002vGRCh38_wm_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.tsv"
#grch38hg002variants = read.table(grch38hg002fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38hg002variantsfiltered <- grch38hg002variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(grch38hg002variantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(grch38hg002variantsfiltered %>% filter(SVTYPE == "DEL"))
#grch38_hg002_length_ofn <- "/home/mkirsche/JasmineT2T/grch38_hg002_hifi.lengths.png"
#plot_length(grch38hg002variants, "GRCh38, HG002", grch38_hg002_length_ofn, T, .92, .85)
#grch38_hg002_length_line_ofn <- "/home/mkirsche/JasmineT2T/grch38_hg002_hifi.lengths_line.png"
#plot_length_line(grch38hg002variants, "GRCh38", grch38_hg002_length_line_ofn, T, .92, .85)


###
### INPUT - CHM13/WINNOWMAP/HIFI
###

#chm13fn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi.specprec.tsv"
#chm13variants = read.table(chm13fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13variantsfiltered = chm13variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(chm13variantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(chm13variantsfiltered %>% filter(SVTYPE == "DEL"))


###
### INDEL BALANCE - CHM13/WINNOWMAP/HIFI
###

#chm13_length_ofn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi.lengths.png"
#plot_length(chm13variants, "CHM13", chm13_length_ofn, T, .92, .85)
#chm13_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi.lengths_line.png"
#plot_length_line(chm13variants, "CHM13", chm13_length_line_ofn, T, .92, .85)

###
### ALLELE FREQUENCIES - CHM13/WINNOWMAP/HIFI
###

#chm13_af_ofn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi.allelefrequency.png"
#plot_allele_frequencies(chm13variants, chm13_af_ofn, "SV Allele Frequencies (CHM13, HiFi, wm)")


###
### INPUT - CHM13/NEWWINNOWMAP/HIFI
###

#chm13newfn <- "/home/mkirsche/JasmineT2T/newwm/chm13_wm_hifi.specprec.tsv"
#chm13newvariants = read.table(chm13newfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13newvariantsfiltered <- chm13newvariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(chm13newvariantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(chm13newvariantsfiltered %>% filter(SVTYPE == "DEL"))


###
### INPUT - CHM13/MINIMAP2/HIFI
###

#chm13mm2fn <- "/home/mkirsche/JasmineT2T/chm13_mm2_hifi.merged.tsv"
#chm13mm2variants = read.table(chm13mm2fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13mm2variantsfiltered <- chm13mm2variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(chm13mm2variantsfiltered %>% filter(SVTYPE == "INS"))
#nrow(chm13mm2variantsfiltered %>% filter(SVTYPE == "DEL"))


###
### INDEL BALANCE - CHM13/MINIMAP2/HIFI
###

#chm13_mm2_length_ofn <- "/home/mkirsche/JasmineT2T/chm13_mm2_hifi.lengths.png"
#plot_length(chm13mm2variants, "CHM13, mm2", chm13_mm2_length_ofn, T, .92, .85)
#chm13_mm2_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_mm2_hifi.lengths_line.png"
#plot_length_line(chm13mm2variants, "CHM13, mm2", chm13_mm2_length_line_ofn, T, .92, .85)


###
### ALLELE FREQUENCIES - CHM13/MINIMAP2/HIFI
###
#chm13_mm2_af_ofn <- "/home/mkirsche/JasmineT2T/chm13_mm2_hifi.allelefrequency.png"
#plot_allele_frequencies(chm13mm2variants, chm13_mm2_af_ofn, "SV Allele Frequencies (CHM13, HiFi, mm2)")


###
### INDEL BALANCE HG002 - CHM13/WINNOWMAP/HIFI
###

#chm13hg002fn <- "/home/mkirsche/JasmineT2T/HG002vCHM13_20200921_wm_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.tsv"
#chm13hg002variants = read.table(chm13hg002fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13hg002variantsfiltered <- chm13hg002variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#nrow(chm13hg002variantsfiltered %>% filter(SVTYPE == "INS" & abs(as.numeric(SVLEN)) >= 50))
#nrow(chm13hg002variantsfiltered %>% filter(SVTYPE == "DEL" & abs(as.numeric(SVLEN)) >= 50))
#chm13_hg002_length_ofn <- "/home/mkirsche/JasmineT2T/chm13_hg002_hifi.lengths.png"
#plot_length(chm13hg002variants, "CHM13, HG002", chm13_hg002_length_ofn, T, .92, .85)
#chm13_hg002_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_hg002_hifi.lengths_line.png"
#plot_length_line(chm13hg002variants, "CHM13", chm13_hg002_length_line_ofn, T, .92, .85)


###
### DISCORDANCE HG002 - GRCH38/WINNOWMAP/HIFI
###

#grch38ashfn <- "/home/mkirsche/JasmineT2T/grch38_wm_hifi_ashtrio.merged.tsv"
#grch38ashvariants = read.table(grch38ashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_ash_wm_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(grch38ashvariants, 'GRCh38', grch38_ash_suppvec_ofn, T, 1700, .9, .85, 'GRCh38 HG002', 20)


###
### DISCORDANCE HG005 - GRCH38/WINNOWMAP/HIFI
###

#grch38hanfn <- "/home/mkirsche/JasmineT2T/grch38_wm_hifi_hantrio.merged.tsv"
#grch38hanvariants = read.table(grch38hanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_han_wm_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(grch38hanvariants, 'GRCh38', grch38_han_suppvec_ofn, T, 1700, .9, .85, 'GRCh38 HG005', 20)


###
### DISCORDANCE HG002 - GRCH38/MINIMAP2/HIFI
###

#grch38mm2ashfn <- "/home/mkirsche/JasmineT2T/grch38_mm2_hifi_ashtrio.merged.tsv"
#grch38mm2ashvariants = read.table(grch38mm2ashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38_mm2_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_ash_mm2_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(grch38mm2ashvariants, 'GRCh38', grch38_mm2_ash_suppvec_ofn, T, 1700, .9, .85, 'GRCh38 HG002, mm2', 20)


###
### DISCORDANCE HG005 - GRCH38/MINIMAP2/HIFI
###

#grch38mm2hanfn <- "/home/mkirsche/JasmineT2T/grch38_mm2_hifi_hantrio.merged.tsv"
#grch38mm2hanvariants = read.table(grch38mm2hanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38_mm2_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_han_mm2_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(grch38mm2hanvariants, 'GRCh38', grch38_mm2_han_suppvec_ofn, T, 1700, .9, .85, 'GRCh38 HG005, mm2', 20)


###
### DISCORDANCE HG002 - GRCH38/MIX/HIFI
###

#grch38intersectashfn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi_ashtrio.specprec.tsv"
#grch38intersectashvariants = read.table(grch38intersectashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38_intersect_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_ash_intersect_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(grch38intersectashvariants, 'GRCh38', grch38_intersect_ash_suppvec_ofn, T, 1700, .9, .85, 'GRCh38 HG002, mm2+wm', 20)


###
### DISCORDANCE HG005 - GRCH38/MIX/HIFI
###

#grch38intersecthanfn <- "/home/mkirsche/JasmineT2T/grch38_intersect_hifi_hantrio.specprec.tsv"
#grch38intersecthanvariants = read.table(grch38intersecthanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#grch38_intersect_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_han_intersect_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(grch38intersecthanvariants, 'GRCh38', grch38_intersect_han_suppvec_ofn, T, 1700, .9, .85, 'GRCh38 HG005, mm2+wm', 20)


###
### DISCORDANCE HG002 - CHM13/MIX/HIFI
###

#chm13intersectashfn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi_ashtrio.specprec.tsv"
#chm13intersectashvariants = read.table(chm13intersectashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_intersect_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ash_intersect_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(chm13intersectashvariants, 'CHM13', chm13_intersect_ash_suppvec_ofn, T, 1700, .9, .85, 'CHM13 HG002, mm2+wm', 20)


###
### DISCORDANCE HG005 - CHM13/MIX/HIFI
###

#chm13intersecthanfn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi_hantrio.specprec.tsv"
#chm13intersecthanvariants = read.table(chm13intersecthanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_intersect_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_han_intersect_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(chm13intersecthanvariants, 'CHM13', chm13_intersect_han_suppvec_ofn, T, 1700, .9, .85, 'CHM13 HG005, mm2+wm', 20)


###
### DISCORDANCE HG002 - GRCH38/MIX/MIX
###

grch38intersectashmixfn <- "/home/mkirsche/JasmineT2T/grch38_ash_fixed.merged.tsv"
grch38intersectashmixvariants = read.table(grch38intersectashmixfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(grch38intersectashmixvariants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1"))
grch38_intersect_ash_mix_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_ash_intersect_mix.suppvec.png"
suppvec_hist_highlightdiscordant(grch38intersectashmixvariants, 'GRCh38', grch38_intersect_ash_mix_suppvec_ofn, T, 2200, .9, .85, 'GRCh38 HG002', 20, 35000)
grch38_intersect_ash_mix_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_ash_intersect_mix.suppvec.svg"
suppvec_hist_highlightdiscordant(grch38intersectashmixvariants, 'GRCh38', grch38_intersect_ash_mix_suppvec_ofn, T, 2200, .88, .85, 'GRCh38 HG002', 20, 35000)


###
### DISCORDANCE HG005 - GRCH38/MIX/MIX
###

grch38intersecthanfn <- "/home/mkirsche/JasmineT2T/grch38_han_fixed.merged.tsv"
grch38intersecthanvariants = read.table(grch38intersecthanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(grch38intersecthanvariants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1" & SUPP_VEC == "100"))
grch38_intersect_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_han_intersect_mix.suppvec.png"
suppvec_hist_highlightdiscordant(grch38intersecthanvariants, 'GRCh38', grch38_intersect_han_suppvec_ofn, T, 2200, .9, .85, 'GRCh38 HG005', 20, 35000)
grch38_intersect_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/grch38_han_intersect_mix.suppvec.svg"
suppvec_hist_highlightdiscordant(grch38intersecthanvariants, 'GRCh38', grch38_intersect_han_suppvec_ofn, T, 2200, .88, .85, 'GRCh38 HG005', 20, 35000)


###
### DISCORDANCE HG002 - CHM13/MIX/MIX
###

chm13intersectashfn <- "/home/mkirsche/JasmineT2T/chm13_ash_fixed.merged.tsv"
chm13intersectashvariants = read.table(chm13intersectashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(chm13intersectashvariants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1"))

chm13_intersect_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ash_intersect_mix.suppvec.png"
suppvec_hist_highlightdiscordant(chm13intersectashvariants, 'CHM13', chm13_intersect_ash_suppvec_ofn, T, 2200, .9, .85, 'CHM13 HG002', 20, 35000, 'HG002', 'HG003', 'HG004')
chm13_intersect_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ash_intersect_mix.suppvec.svg"
suppvec_hist_highlightdiscordant(chm13intersectashvariants, 'CHM13', chm13_intersect_ash_suppvec_ofn, T, 2200, .88, .85, 'CHM13 HG002', 20, 35000, 'HG002', 'HG003', 'HG004')


###
### DISCORDANCE HG005 - CHM13/MIX/MIX
###

chm13intersecthanfn <- "/home/mkirsche/JasmineT2T/chm13_han_fixed.merged.tsv"
chm13intersecthanvariants = read.table(chm13intersecthanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
nrow(chm13intersecthanvariants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1"))
chm13_intersect_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_han_intersect_mix.suppvec.png"
suppvec_hist_highlightdiscordant(chm13intersecthanvariants, 'CHM13', chm13_intersect_han_suppvec_ofn, T, 2200, .9, .85, 'CHM13 HG005', 20, 35000, 'HG005', 'HG006', 'HG007')
chm13_intersect_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_han_intersect_mix.suppvec.svg"
suppvec_hist_highlightdiscordant(chm13intersecthanvariants, 'CHM13', chm13_intersect_han_suppvec_ofn, T, 2200, .88, .85, 'CHM13 HG005', 20, 35000, 'HG005', 'HG006', 'HG007')


###
### DISCORDANCE HG002 BY GENOMIC CONTEXT - CHM13/WINNOWMAP/HIFI
###

#chm13ashannofn <- "/home/mkirsche/JasmineT2T/chm13_ash_anno_annotated.tsv"
#chm13ashannovariants = read.table(chm13ashannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_ashanno_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ashanno_wm_hifi.suppvec.png"
#filtered <- chm13ashannovariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#colnames(filtered)
#div <- 40
#all <- suppvec_hist_highlightdiscordant2(filtered, "All Variants", div)
#exon <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_EXON == 1), "Exon", div)
#gene <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_GENE == 1), "Gene", div)
#synteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_SYNTENY_CHM13 == 1), "Syntenic", div)
#nonsynteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), "Non-Syntenic", div + 10)
#ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
#ggsave(chm13_ashanno_suppvec_ofn, width=20, height=16)


###
### DISCORDANCE HG002 BY GENOMIC CONTEXT - CHM13/MIX/HIFI
###

chm13intersectashannofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_ashtrio_anno_annotated.tsv"
chm13intersectashannovariants = read.table(chm13intersectashannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
chm13_intersect_ashanno_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ashanno_intersect.suppvec.png"
filtered <- chm13intersectashannovariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
colnames(filtered)
div <- 40
all <- suppvec_hist_highlightdiscordant2(filtered, "All Variants", div)
exon <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_EXON == 1), "Exon", div)
gene <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_GENE == 1), "Gene", div)
synteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_SYNTENY_CHM13 == 1), "Syntenic", div)
nonsynteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), "Non-Syntenic", div + 10)
ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
ggsave(chm13_intersect_ashanno_suppvec_ofn, width=20, height=16)
chm13_intersect_ashanno_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ashanno_intersect.suppvec.svg"
ggsave(chm13_intersect_ashanno_suppvec_ofn, width=20, height=16)


###
### DISCORDANCE HG002 - CHM13/WINNOWMAP/HIFI
###

#chm13ashfn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi_ashtrio.merged.tsv"
#chm13ashvariants = read.table(chm13ashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ash_wm_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(chm13ashvariants, 'CHM13', chm13_ash_suppvec_ofn, T, 1700, .9, .85, 'CHM13 HG002', 20)


###
### DISCORDANCE HG005 BY GENOMIC CONTEXT - CHM13/WINNOWMAP/HIFI
###

#chm13hanannofn <- "/home/mkirsche/JasmineT2T/chm13_han_anno_annotated.tsv"
#chm13hanannovariants = read.table(chm13hanannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_hananno_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_hananno_wm_hifi.suppvec.png"
#filtered <- chm13hanannovariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#colnames(filtered)
#div <- 40
#all <- suppvec_hist_highlightdiscordant2(filtered, "All Variants", div)
#exon <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_EXON == 1), "Exon", div)
#gene <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_GENE == 1), "Gene", div)
#synteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_SYNTENY_CHM13 == 1), "Syntenic", div)
#nonsynteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), "Non-Syntenic", div + 10)
#ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
#ggsave(chm13_hananno_suppvec_ofn, width=20, height=16)


###
### DISCORDANCE HG005 BY GENOMIC CONTEXT - CHM13/MIX/HIFI
###

chm13intersecthanannofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hantrio_anno_annotated.tsv"
chm13intersecthanannovariants = read.table(chm13intersecthanannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
chm13_intersect_hananno_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_hananno_intersect.suppvec.png"
filtered <- chm13intersecthanannovariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
colnames(filtered)
div <- 40
all <- suppvec_hist_highlightdiscordant2(filtered, "All Variants", div)
exon <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_EXON == 1), "Exon", div)
gene <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_GENE == 1), "Gene", div)
synteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_SYNTENY_CHM13 == 1), "Syntenic", div)
nonsynteny <- suppvec_hist_highlightdiscordant2(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), "Non-Syntenic", div + 10)
ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
ggsave(chm13_intersect_hananno_suppvec_ofn, width=20, height=16)
chm13_intersect_hananno_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_hananno_intersect.suppvec.svg"
ggsave(chm13_intersect_hananno_suppvec_ofn, width=20, height=16)


###
### DISCORDANCE HG005 - CHM13/WINNOWMAP/HIFI
###

#chm13hanfn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi_hantrio.merged.tsv"
#chm13hanvariants = read.table(chm13hanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_han_wm_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(chm13hanvariants, 'CHM13', chm13_han_suppvec_ofn, T, 1700, .9, .85, 'CHM13 HG005', 20)


###
### DISCORDANCE HG002 - CHM13/MINIMAP2/HIFI
###

#chm13mm2ashfn <- "/home/mkirsche/JasmineT2T/chm13_mm2_hifi_ashtrio.merged.tsv"
#chm13mm2ashvariants = read.table(chm13mm2ashfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_mm2_ash_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_ash_mm2_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(chm13mm2ashvariants, 'CHM13', chm13_mm2_ash_suppvec_ofn, T, 1700, .9, .85, 'CHM13 HG002, mm2', 20)

###
### DISCORDANCE HG005 - CHM13/MINIMAP2/HIFI
###

#chm13mm2hanfn <- "/home/mkirsche/JasmineT2T/chm13_mm2_hifi_hantrio.merged.tsv"
#chm13mm2hanvariants = read.table(chm13mm2hanfn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_mm2_han_suppvec_ofn <- "/home/mkirsche/JasmineT2T/chm13_han_mm2_hifi.suppvec.png"
#suppvec_hist_highlightdiscordant(chm13mm2hanvariants, 'CHM13', chm13_mm2_han_suppvec_ofn, T, 1700, .9, .85, 'CHM13 HG005, mm2', 20)


###
### SIZE BY GENOMIC CONTEXT - CHM13/WINNOWMAP/HIFI
###

#chm13annofn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi_anno_annotated.tsv"
#chm13annovariants = read.table(chm13annofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
#chm13_anno_types_ofn <- "/home/mkirsche/JasmineT2T/chm13_anno_wm_hifi.types.png"
#filtered <- chm13annovariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#colnames(filtered)
#all <- plot_table(filtered, "All Variants")
#exon <- plot_table(filtered %>% filter(INTERSECTS_EXON == 1), "Exon")
#gene <- plot_table(filtered %>% filter(INTERSECTS_GENE == 1), "Gene")
#synteny <- plot_table(filtered %>% filter(INTERSECTS_SYNTENY_CHM13 == 1), "Syntenic")
#nonsynteny <- plot_table(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), "Non-Syntenic")
#ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
#ggsave(chm13_anno_types_ofn, width=20, height=16)
#chm13_anno_types_ofn <- "/home/mkirsche/JasmineT2T/chm13_anno_wm_hifi.types.svg"
#ggsave(chm13_anno_types_ofn, width=20, height=16)


#highlighted_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi_anno.lengths_line.png"
#plot_length_line_highlight(chm13annovariants, "CHM13", highlighted_length_line_ofn, T, .92, .85)
#highlighted_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_wm_hifi_anno.lengths_line.svg"
#plot_length_line_highlight(chm13annovariants, "CHM13", highlighted_length_line_ofn, T, .92, .85)

###
### SIZE BY GENOMIC CONTEXT - CHM13/MIX/HIFI
###

chm13intersectanno2fn <- "/home/mkirsche/JasmineT2T/chm13_intersect_anno2_annotated.tsv"
chm13intersectanno2variants = read.table(chm13intersectanno2fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
filtered2 <- chm13intersectanno2variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
nrow(filtered2 %>% filter(INTERSECTS_NOVEL_CHM13 == 1 & INTERSECTS_NONSYNTENY_CHM13 == 1 & SVTYPE== "DEL"))

chm13intersectannofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_anno_annotated.tsv"
chm13intersectannovariants = read.table(chm13intersectannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
chm13_intersect_anno_types_ofn <- "/home/mkirsche/JasmineT2T/chm13_anno_intersect.types.png"
chm13_intersect_anno_repeat_types_ofn <- "/home/mkirsche/JasmineT2T/chm13_anno_repeat_intersect.types.png"
filtered <- chm13intersectannovariants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)

chm13_intersect_af_ns_ofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_hifi_nonsyntenic.allelefrequency.png"
plot_allele_frequencies(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), chm13_intersect_af_ns_ofn, "SV Allele Frequencies (CHM13 Non-syntenic)", 12000)
nrow(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1))
nrow(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1 & SVTYPE == "INS"))
nrow(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1 & SVTYPE == "DEL"))

colnames(filtered)
all <- plot_table(filtered, "All Variants")
exon <- plot_table(filtered %>% filter(INTERSECTS_EXON == 1), "Exon")
gene <- plot_table(filtered %>% filter(INTERSECTS_GENE == 1), "Gene")
synteny <- plot_table(filtered %>% filter(INTERSECTS_SYNTENY_CHM13 == 1), "Syntenic")
nonsynteny <- plot_table(filtered %>% filter(INTERSECTS_NONSYNTENY_CHM13 == 1), "Non-Syntenic")
line <- plot_table(filtered %>% filter(INTERSECTS_LINE == 1), "LINE Repeat")
sine <- plot_table(filtered %>% filter(INTERSECTS_SINE == 1), "SINE Repeat")
satellite <- plot_table(filtered %>% filter(INTERSECTS_SATELLITE == 1), "Satellite Repeat")
lowcomplexity <- plot_table(filtered %>% filter(INTERSECTS_LOW_COMPLEXITY == 1), "Low Complexity Repeat")
simple <- plot_table(filtered %>% filter(INTERSECTS_SIMPLE_REPEAT == 1), "Simple Repeat")
centromere <- plot_table(filtered %>% filter(INTERSECTS_CENTROMERE == 1), "Centromere")

ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
ggsave(chm13_intersect_anno_types_ofn, width=20, height=16)

ggarrange(line, sine, satellite, lowcomplexity, simple, centromere,nrow = 2, ncol=3)
ggsave(chm13_intersect_anno_repeat_types_ofn, width=20, height=16)

chm13_intersect_anno_types_ofn <- "/home/mkirsche/JasmineT2T/chm13_anno_intersect.types.svg"
chm13_intersect_anno_repeat_types_ofn <- "/home/mkirsche/JasmineT2T/chm13_anno_repeat_intersect.types.svg"

ggarrange(all,exon, gene, synteny, nonsynteny,nrow = 2, ncol=3)
ggsave(chm13_intersect_anno_types_ofn, width=20, height=16)

ggarrange(line, sine, satellite, lowcomplexity, simple, centromere,nrow = 2, ncol=3)
ggsave(chm13_intersect_anno_repeat_types_ofn, width=20, height=16)

highlighted_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line.png"
plot_length_line_highlight(chm13intersectannovariants, "CHM13", highlighted_length_line_ofn, T, .92, .85)
highlighted_length_line_ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line.svg"
plot_length_line_highlight(chm13intersectannovariants, "CHM13", highlighted_length_line_ofn, T, .92, .85)

### 
### GENOMIC POSITIONS - CHM13/WINNOWMAP/HIFI
###

goodlist <- paste('chr', c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', 
                           '21', '22', "X", "Y"), sep = "")
binsize <- 100000
chm13variants$STARTBIN <- floor(as.numeric(chm13variants$POS) / binsize) * binsize
groups <- chm13variants %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
unique(groups$CHROM)
groups$CHROM = factor(groups$CHROM, levels=goodlist)
groups$logn <- log2(groups$n + 1)
outfile <- "/home/mkirsche/JasmineT2T/chm13_chr_svs.png"
ggplot(groups %>% filter(CHROM %in% goodlist), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity") +
  xlab("Genomic Position") +
  ylab("SVs per 100kbp bin (log2)") +
  ggtitle("SV Density (CHM13)") +
  facet_wrap(~ CHROM, ncol = 5, scales = "free") +
  theme(plot.title = element_text(hjust = .5)) 
ggsave(outfile, width = 12, height = 12)

### 
### HG002 GENOMIC POSITIONS - CHM13/WINNOWMAP/HIFI
###

#chm13hg002variants <- chm13hg002variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
#binsize <- 100000
#chm13hg002variants$STARTBIN <- floor(as.numeric(chm13hg002variants$POS) / binsize) * binsize
#groups <- chm13hg002variants %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
#unique(groups$CHROM)
#groups$CHROM = factor(groups$CHROM, levels=goodlist)
#groups$logn <- log2(groups$n + 1)
#outfile <- "/home/mkirsche/JasmineT2T/chm13_hg002_chr_svs.png"
#ggplot(groups %>% filter(CHROM %in% goodlist), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity") +
#  xlab("Genomic Position") +
#  ylab("SVs per 100kbp bin (log2)") +
#  ggtitle("SV Density (CHM13)") +
#  facet_wrap(~ CHROM, ncol = 5, scales = "free") +
#  theme(plot.title = element_text(hjust = .5)) 
#ggsave(outfile, width = 12, height = 12)


### 
### GENOMIC POSITIONS - GRCH38/WINNOWMAP/HIFI
###

grch38variants$STARTBIN <- floor(as.numeric(grch38variants$POS) / binsize) * binsize
groups <- grch38variants %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
unique(groups$CHROM)
groups$CHROM = factor(groups$CHROM, levels=goodlist)
groups$logn <- log2(groups$n + 1)
outfile <- "/home/mkirsche/JasmineT2T/grch38_chr_svs.png"
ggplot(groups %>% filter(CHROM %in% goodlist), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity") +
  xlab("Genomic Position") +
  ylab("SVs per 100kbp bin (log2)") +
  ggtitle("SV Density (GRCh38)") +
  facet_wrap(~ CHROM, ncol = 5, scales = "free") +
  theme(plot.title = element_text(hjust = .5)) 
ggsave(outfile, width = 12, height = 12)

### 
### HG002 GENOMIC POSITIONS - GRCH38/WINNOWMAP/HIFI
###

grch38hg002variants <- grch38hg002variants %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
grch38hg002variants$STARTBIN <- floor(as.numeric(grch38hg002variants$POS) / binsize) * binsize
groups <- grch38hg002variants %>% group_by(CHROM, STARTBIN, SVTYPE) %>% tally()
groups
unique(groups$CHROM)
groups$CHROM = factor(groups$CHROM, levels=goodlist)
groups$logn <- log2(groups$n + 1)
outfile <- "/home/mkirsche/JasmineT2T/grch38_hg002_chr_svs.png"
ggplot(groups %>% filter(CHROM %in% goodlist), aes(x = STARTBIN, y = logn, fill = SVTYPE)) + geom_bar( stat = "identity") +
  xlab("Genomic Position") +
  ylab("SVs per 100kbp bin (log2)") +
  ggtitle("SV Density (GRCh38)") +
  facet_wrap(~ CHROM, ncol = 5, scales = "free") +
  theme(plot.title = element_text(hjust = .5)) 
ggsave(outfile, width = 12, height = 12)


###
### PAIRWISE CLUSTERING
###

library(pheatmap)
library(grid)
library(RColorBrewer)

get_pairwise <- function(fn) {
  variants = read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
  variants <- variants %>% filter(IS_SPECIFIC == "1" & IS_PRECISE == "1")
  nrow(variants)
  length <- max(nchar(variants$SUPP_VEC))
  mat <- array(dim = c(length, length))
  for(i in 1:length) {
    for(j in 1:length) {
      #mat[i, j] = nrow(variants %>% filter(charAt(as.character(filter$SUPP_VEC, i)) == '1'))
      mat[i, j] = sum(substr(variants$SUPP_VEC,i,i) == "1" & substr(variants$SUPP_VEC,j,j) == "1")
      if(i == j) {
        #mat[i, j] = 0
      }
    } 
  }
  return(mat)
}
mat <- get_pairwise(chm13fn)

# Read in sample info
fn <- "/home/mkirsche/JasmineT2T/filelist.txt"
samples <- read.table(fn, sep = "\t", header = T, colClasses = 'character')
samples$SAMPLE
samples$Label <- samples$SAMPLE
samples
rownames(mat) <- samples$Label
colnames(mat) <- samples$Label
#rownames(mat_clr100) <- samples$Label
#colnames(mat_clr100) <- samples$Label

samples$COVERAGE_NUM <- as.numeric(samples$COVERAGE)
samples$NUMVARS <- as.numeric(samples$NUMVARS)

mat
mat_colors <- list(Superpopulation = brewer.pal(5, "Set1"))
names(mat_colors$Superpopulation) <- unique(samples$SUPERPOPULATION)
mat_colors
mat_col <- data.frame(Superpopulation = samples$SUPERPOPULATION)
rownames(mat_col) <- rownames(mat)
mat_col
outfile <- "/home/mkirsche/JasmineT2T/chm13pairwise.png"
png(file=outfile, width = 1024, height = 1024)
pheatmap(mat,show_rownames=TRUE,show_colnames=TRUE, annotation_colors = mat_colors, annotation_col = mat_col, annotation_row = mat_col, drop_levels = TRUE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean", fontsize = 16)
dev.off()

outfile <- "/home/mkirsche/JasmineT2T/chm13pairwise.svg"
#svg(file=outfile, width = 1024, height = 1024)
Cairo(1024, 1024, file = outfile, type = "svg", bg = "white", dpi = 80)
pheatmap(mat,show_rownames=TRUE,show_colnames=TRUE, annotation_colors = mat_colors, annotation_col = mat_col, annotation_row = mat_col, drop_levels = TRUE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean", fontsize = 12)
dev.off()

grch38mat <- get_pairwise(grch38fn)
grch38mat
rownames(grch38mat) <- samples$Label
colnames(grch38mat) <- samples$Label
mat_colors <- list(Superpopulation = brewer.pal(5, "Set1"))
names(mat_colors$Superpopulation) <- unique(samples$SUPERPOPULATION)
mat_colors
mat_col <- data.frame(Superpopulation = samples$SUPERPOPULATION)
rownames(mat_col) <- rownames(grch38mat)
mat_col
outfile <- "/home/mkirsche/JasmineT2T/grch38pairwise.png"
png(file=outfile, width = 1024, height = 1024)
pheatmap(grch38mat,show_rownames=TRUE,show_colnames=TRUE, annotation_colors = mat_colors, annotation_col = mat_col, annotation_row = mat_col, drop_levels = TRUE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean", fontsize = 16)
dev.off()

outfile <- "/home/mkirsche/JasmineT2T/grch38pairwise.svg"
#svg(file=outfile, width = 1024, height = 1024)
Cairo(1024, 1024, file = outfile, type = "svg", bg = "white", dpi = 80)
pheatmap(grch38mat,show_rownames=TRUE,show_colnames=TRUE, annotation_colors = mat_colors, annotation_col = mat_col, annotation_row = mat_col, drop_levels = TRUE,
         scale = "none",clustering_method="ward.D2",
         clustering_distance_cols="euclidean", fontsize = 12)
dev.off()
