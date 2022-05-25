library(ggplot2)
library(ggpattern)
library(dplyr)

chm13intersectannofn <- "/home/mkirsche/JasmineT2T/chm13_intersect_anno_annotated.tsv"
df = read.table(chm13intersectannofn, comment.char = "#", sep = "\t", header = T, stringsAsFactors=FALSE, colClasses = 'character')
ofn <- "/home/mkirsche/JasmineT2T/chm13_mix_hifi_anno.lengths_line_pattern.png"
plot_length_line_highlight(chm13intersectannovariants, "CHM13", highlighted_length_line_ofn, T, .92, .85)

df <- df %>% filter(IS_SPECIFIC == 1 & IS_PRECISE == 1)
  summarized <- df %>% group_by(SVTYPE, INTERSECTS_NONSYNTENY_CHM13) %>% summarise(counts=n()) %>% filter(SVTYPE == "INS" | SVTYPE == "DEL")
  summarized$CATEGORY <- paste(summarized$SVTYPE, "_", summarized$INTERSECTS_NONSYNTENY_CHM13, sep = "")
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "INS_0", "Syntenic Insertion",summarized$CATEGORY)
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "INS_1", "Non-syntenic Insertion",summarized$CATEGORY)
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "DEL_0", "Syntenic Deletion",summarized$CATEGORY)
  summarized$CATEGORY <- ifelse(summarized$CATEGORY == "DEL_1", "Non-syntenic Deletion",summarized$CATEGORY)
  summarized
  ggplot(data = summarized %>% filter(SVTYPE == "INS" | SVTYPE == "DEL"), aes(x = SVTYPE, y = counts, fill = SVTYPE, pattern = INTERSECTS_NONSYNTENY_CHM13)) +
    geom_bar_pattern(position = position_dodge(preserve = "single"), stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) + 
    scale_fill_manual(values = c('#CC6677', '#44AA99')) +
    scale_pattern_manual(values = c(1 = "stripe", 0 = "none")) +
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
    ylim(0, 70000) +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))
  ggsave(ofn, width= 3, height = 8)


         
         
      
