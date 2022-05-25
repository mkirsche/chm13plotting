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

indir <- "/home/mkirsche/JasmineT2T/ReadLength/"
wildcard <- "*ONT*.txt"
files <- Sys.glob(paste(indir, wildcard, sep = ""))


names = c("Length", "SAMPLE", "TECH")
all_data <- data.frame(setNames(rep(list(NA), length(names)), names))
all_data$File = c()
all_data$Total = c()
all_data$Discordant = c()
md_file <- "/home/mkirsche/JasmineT2T/ReadLength/HG002.ONT.readslength.txt"
for(md_file in files) {
  md_file
  md_df = read.table(md_file, comment.char = "#", sep = "\t", header = F)
  colnames(md_df) <- c("Length")
  min(md_df$Length)
  md_df <- md_df %>% filter(as.numeric(Length) >= 1000)
  min(md_df$Length)
  
  samplename <- gsub(".ONT.readslength.txt", '', substring(md_file, str_length(indir) + 1))
  md_df$SAMPLE <- samplename
  md_df$TECH <- "ONT"
  #nrow(md_df)
  
  #md_df_filter <- md_df %>% filter(SUPP_VEC == "100")
  
  #disccount <- nrow(md_df_filter)
  #disccount
  #totalcount <- nrow(md_df)
  
  #cur <- data.frame(md_file, totalcount, disccount)
  #cur
  #names(cur) <- c("File", "Total", "Discordant")
  #md_data <- rbind(md_data, cur)
  ggplot(md_df, aes(x = log2(as.numeric(Length)))) + geom_density()
  ggsave(paste(indir, samplename, '_', "ONT", "_hist.png", sep = ""), width = 8, height = 8)
  
  all_data <- rbind(all_data, sample_n(md_df, nrow(md_df)/20))
  
}
all_data <- all_data %>% filter(!is.na(SAMPLE))# & as.numeric(Length) >= 1000)
min(all_data$Length)

#ggplot(all_data, aes(x = log2(as.numeric(Length)), color = SAMPLE, fill = SAMPLE)) + geom_density(alpha = .1)
#ggsave(paste(indir, "all", "_ONT_hist.pdf", sep = ""), device = "pdf", width = 8, height = 8)
ggplot(all_data %>% filter(!is.na(SAMPLE)), aes(x = (as.numeric(Length)), y = SAMPLE, color = SAMPLE, fill = SAMPLE)) + geom_violin() + xlab("Length") + ylab("Sample") + theme(legend.position = "none") +
  stat_summary(fun.data=mean_sdl, geom="pointrange", color = "black", fun.args = list(mult = 1)) + xlim(-10000, 100000)
ggsave(paste(indir, "all", "_ONT_violin.svg", sep = ""), width = 8, height = 8)
ggsave(paste(indir, "all", "_ONT_violin.png", sep = ""), width = 8, height = 8)

  summ_df <- all_data %>% filter(!is.na(SAMPLE)) %>%
  group_by(SAMPLE) %>%
  summarize(mean=mean(Length))
summ_df$mean
mean(summ_df$mean)
tmp <- all_data
###
### HiFi
###
wildcard <- "*HiFi*.txt"
files <- Sys.glob(paste(indir, wildcard, sep = ""))
names = c("Length", "SAMPLE", "TECH")
all_data <- data.frame(setNames(rep(list(NA), length(names)), names))
all_data$File = c()
all_data$Total = c()
all_data$Discordant = c()
for(md_file in files) {
  md_file
  md_df = read.table(md_file, comment.char = "#", sep = "\t", header = F)
  colnames(md_df) <- c("Length")
  samplename <- gsub(".HiFi.readslength.txt", '', substring(md_file, str_length(indir) + 1))
  md_df$SAMPLE <- samplename
  md_df$TECH <- "HiFi"
  #nrow(md_df)
  
  #md_df_filter <- md_df %>% filter(SUPP_VEC == "100")
  
  #disccount <- nrow(md_df_filter)
  #disccount
  #totalcount <- nrow(md_df)
  
  #cur <- data.frame(md_file, totalcount, disccount)
  #cur
  #names(cur) <- c("File", "Total", "Discordant")
  #md_data <- rbind(md_data, cur)
  ggplot(md_df, aes(x = log2(as.numeric(Length)))) + geom_density()
  ggsave(paste(indir, samplename, '_', "HiFi", "_violin.png", sep = ""), width = 8, height = 8)
  
  all_data <- rbind(all_data, sample_n(md_df, nrow(md_df)/20))
  
}

tmp <- rbind(tmp, all_data)  

summ_df <- all_data %>%
  group_by(SAMPLE) %>%
  summarize(mean=mean(log2(Length)))
ggplot(all_data, aes(x = as.numeric(Length), color = SAMPLE, fill = SAMPLE)) + geom_density(alpha = .1)
ggsave(paste(indir, "all", "_HiFi_hist.pdf", sep = ""), device = "pdf", width = 8, height = 8)

ggplot(all_data %>% filter(!is.na(SAMPLE)), aes(x = as.numeric(Length), y = SAMPLE, color = SAMPLE, fill = SAMPLE)) + geom_violin() + xlab("Length") + ylab("Sample") + theme(legend.position = "none") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color = "black")
ggsave(paste(indir, "all", "_HiFi_violin.svg", sep = ""), width = 8, height = 8)
ggsave(paste(indir, "all", "_HiFi_violin.png", sep = ""), width = 8, height = 8)
summ_df <- all_data %>% filter(!is.na(SAMPLE)) %>%
  group_by(SAMPLE) %>%
  summarize(mean=mean(Length))

mean(summ_df$mean)

ggplot(tmp %>% filter(!is.na(SAMPLE)), aes(x = as.numeric(Length), y = SAMPLE, color = SAMPLE, fill = SAMPLE)) + geom_violin() + xlab("Length") + ylab("Sample") + theme(legend.position = "none") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", color = "black") + facet_grid(cols=vars(TECH), scales = "free_x")
ggsave(paste(indir, "all", "_BothTechs_violin.svg", sep = ""), width = 16, height = 8)

