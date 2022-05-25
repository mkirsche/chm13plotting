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

fn <- "/home/mkirsche/JasmineT2T/abcov_summary.bed"
df <- read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
df
colnames(df)
df$Mean <- as.numeric(df$Mean)
df$StandardDeviation <- as.numeric(df$StandardDeviation)
df$Aligner <-  ifelse(df$Aligner == "mm2", "minimap2", df$Aligner)
df$Aligner <-  ifelse(df$Aligner == "wm", "Winnowmap", df$Aligner)
summary <- df %>%
  group_by(Reference, Technology, Aligner) %>%
  summarise_at(vars(Mean, StandardDeviation), funs(mean(.,rm.NA = T)))
summary
colnames(df)

colorpalette <- c(CHM13 = "#88CCEE", GRCh38 = "#AA4499")

df$Reference <- factor(df$Reference, levels = c("GRCh38", "CHM13"))

ggplot(data = df, aes(x = Reference, y = StandardDeviation, fill = Reference)) + #geom_bar(aes(x = Reference, y = Mean, fill = Reference),stat = "identity") +
  geom_boxplot() +
  ylab("Abnormal Coverage Standard Deviation") +
    geom_jitter(color="black", size=0.4) +
    #theme_ipsum() +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.position = "none",
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
    ) +
  scale_fill_manual(values = colorpalette) +
    xlab("") +
  facet_grid(cols = vars(Aligner), rows = vars(Technology), scales = "free_y")
ggsave("/home/mkirsche/JasmineT2T/abcov.png", width = 9,height = 8)
ggsave("/home/mkirsche/JasmineT2T/abcov.svg", width = 9,height = 8)


fn <- "/home/mkirsche/JasmineT2T/abcov_all_summary.bed"
df <- read.table(fn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
df
colnames(df)
df$Mean <- as.numeric(df$Mean)
df$StandardDeviation <- as.numeric(df$StandardDeviation)
df$Aligner <-  ifelse(df$Aligner == "mm2", "minimap2", df$Aligner)
df$Aligner <-  ifelse(df$Aligner == "wm", "Winnowmap", df$Aligner)
summary <- df %>%
  group_by(Reference, Technology, Aligner) %>%
  summarise_at(vars(Mean, StandardDeviation), funs(mean(.,rm.NA = T)))
summary
colnames(df)

colorpalette <- c(CHM13 = "#88CCEE", GRCh38 = "#AA4499")
df$Reference <- factor(df$Reference, levels = c("GRCh38", "CHM13"))
ggplot(data = df, aes(x = Reference, y = StandardDeviation, fill = Reference)) + #geom_bar(aes(x = Reference, y = Mean, fill = Reference),stat = "identity") +
  geom_boxplot() +
  ylab("Coverage Standard Deviation") +
  geom_jitter(color="black", size=0.4) +
  #theme_ipsum() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.position = "none",
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
  ) +
  scale_fill_manual(values = colorpalette) +
  xlab("") +
  facet_grid(cols = vars(Aligner), rows = vars(Technology), scales = "free_y")
ggsave("/home/mkirsche/JasmineT2T/abcovall.png", width = 9,height = 8)
ggsave("/home/mkirsche/JasmineT2T/abcovall.svg", width = 9,height = 8)


fn <- "/home/mkirsche/JasmineT2T/LR.satools.SN.stats.csv"
df <- read.table(fn, comment.char = "#", sep = ",", header = T, stringsAsFactors=FALSE, colClasses = 'character')
df$tech <-  ifelse(df$tech == "PBCCS", "HiFi", df$tech)
df$ref <-  ifelse(df$ref == "GRCh38-noAD", "GRCh38", df$ref)
df$ref <-  ifelse(df$ref == "CHM13_20200921", "CHM13", df$ref)
df$aligner <-  ifelse(df$aligner == "mm2", "minimap2", df$aligner)
df$aligner <-  ifelse(df$aligner == "wm", "Winnowmap", df$aligner)
df$ref <- factor(df$ref, levels = c("GRCh38", "CHM13"))

colnames(df)
ggplot(df, aes(x=ref, y=as.numeric(error.rate), fill=ref)) +
  geom_boxplot() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4) +
  #theme_ipsum() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.position = "none",
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
  ) +
  ylab("Mapping Error Rate") +
  scale_fill_manual(values = colorpalette) +
  xlab("") +
  facet_grid(cols = vars(aligner), rows = vars(tech), scales = "free_y")
ggsave("/home/mkirsche/JasmineT2T/errorrate.png", width = 9, height = 8)
ggsave("/home/mkirsche/JasmineT2T/errorrate.svg", width = 9, height = 8)

