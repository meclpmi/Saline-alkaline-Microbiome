# Project: Saline‑Alkaline Genome Collection (SAGC)
# Figure: Fig1 b‑h (Panel a manually drawn in Illustrator)
# R version: 4.2.2
# Packages: tidyverse, sf, patchwork
# Input: data/SGB_info.csv
# Output: result/*.jpg, result/*.pdf, result/*.tiff
setwd("F:/Saline-alkaline-Microbiome/Metagenome/MAG_stat")
library(tidyverse)
library(sf)
library(patchwork)
library(ggstar)
detach("package:plyr", unload = TRUE)
library(dplyr)
df <- read.csv("data/SGB_info.csv", header = TRUE, fileEncoding = "GBK")
df$Quality.level <- factor(df$Quality.level, levels = c("High-quality","Medium-quality"))
df_agg <- df %>%
  group_by(Lon, Lat, Quality.level) %>%
  summarise(count = n(), .groups = "drop")
df_agg <- df_agg[!is.na(df_agg$Lon) & !is.na(df_agg$Lat), ]
library(ggplot2)
library(rnaturalearth)

world_sf <- ne_countries(scale = 50, returnclass = "sf")
p <- ggplot() +
  geom_sf(data = world_sf,
          fill = "gray80",    # 大陆浅灰
          color = "gray80",   # 海岸线
          linewidth = 0.1) +
  theme_void()
p

p_b <- p +
  geom_point(
    data = df_agg,
    aes(
      x = Lon,
      y = Lat,
      color = Quality.level,
      size = count
    ),
    shape = 1,        # shape=1：纯空心圆圈（只有外圈，内部空白，和目标图完全一致）
    fill = NA,
    stroke = 1,
    position = position_jitter(width = 3, height = 3) # 极小抖动，防止完全重叠，不漂移地理位置
  ) +
  scale_color_manual(
    values = c("Medium-quality" = "#FDE725", "High-quality" = "#5DC863"),
    name = "Quality Level",
    na.value = "grey50"
  ) +
  scale_size_continuous(
    range = c(1, 8),    # 圆圈尺寸范围，如果圈太大就把8改小到5‑6
    name = "Count",
    breaks = c(1,10,100,200),
    labels = c("1","10","100","200")
  ) +
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE, order = 1),
    size  = guide_legend(nrow = 2, byrow = TRUE, order = 2)
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1, "cm"),
    panel.border = element_blank()
  )
p_b
ggsave("result/Fig1b_map.tiff", p_b, width = 8, height = 5, dpi = 300)
ggsave("result/Fig1b_map.jpg", p_b, width = 8, height = 5, dpi = 300)
ggsave("result/Fig1b_map.pdf", p_b, width = 8, height = 5, dpi = 300)

library(ggplot2)
library(dplyr)
library(stringr)

df <- read.csv("data/SGB_info.csv", header = TRUE, fileEncoding = "GBK")

df$cluster <- str_extract(df[,1], "^[^.]+")
df_count <- df %>%
  group_by(Database.type, cluster, Quality.level) %>%
  summarise(MAG_count = n(), .groups = "drop")

df_count$Database <- factor(df_count$Database.type,
                            levels = c("The Aral Sea Basin", "Public"))
df_count$Quality <- factor(df_count$Quality.level,
                           levels = c("Medium-quality", "High-quality"))

ggplot(df_count, aes(
  y = Database,
  x = MAG_count,
  fill = Quality,
  group = interaction(Database, Quality)
)) +
  geom_violin(
    position = position_dodge(width = 0.75),
    scale = "width", trim = TRUE,
    color = NA, alpha = 0.8
  ) +
  geom_jitter(
    color = "black", size = 0.6, alpha = 0.6,
    position = position_jitterdodge(
      jitter.width = 0,
      jitter.height = 0.16,
      dodge.width = 0.75
    )
  ) +
  geom_boxplot(
    width = 0.13, outlier.shape = NA, fill = "white",
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(
    values = c("Medium-quality"="#FDE725", "High-quality"="#5DC863"),
    name = "Quality Level"
  ) +
  labs(x = "Number of MAGs per SGB", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    panel.border = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  )
ggsave("result/Fig1b_volin.tiff", width = 6, height = 5, dpi = 300)
ggsave("result/Fig1b_volin.jpg", width = 6, height = 5, dpi = 300)
ggsave("result/Fig1b_volin.pdf", width = 6, height = 5, dpi = 300)

library(tidyverse)
library(ggridges)

df_raw <- read.csv("data/SGB_info.csv", header = TRUE, fileEncoding = "GBK")

df_long <- df_raw %>%
  pivot_longer(
    cols = c(completeness, contamination, Genome.size, N50, strain_heterogeneity),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  # 固定指标上下顺序（和原图从上到下保持一致！）
  mutate(
    Metric = factor(Metric,
                    levels = c("completeness",
                               "contamination",
                               "Genome.size",
                               "N50",
                               "strain_heterogeneity"))
  )

p <- ggplot(df_long, aes(x = Value, y = Quality.level, fill = Quality.level)) +
  geom_density_ridges(
    scale = 0.7,
    alpha = 0.7,
    color = NA,
    rel_min_height = 0.01
  ) +
  geom_jitter(
    aes(color = Quality.level),
    size = 1.2,
    alpha = 1,
    width = 0,
    height = 0.15
  ) +
  geom_boxplot(
    width = 0.12,
    fill = NA,
    color = "black",
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = c(
    "High-quality" = "#5DC863",
    "Medium-quality" = "#FDE725"
  )) +
  scale_color_manual(values = c(
    "High-quality" = "#5DC863",
    "Medium-quality" = "#FDE725"
  )) +
  facet_wrap(~Metric, ncol = 1, scales = "free_x", strip.position = "bottom") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0.8, "cm"),
    panel.border = element_blank()
  )
p
ggsave("result/Fig1c.tiff", p, width = 7, height = 12, dpi = 600)
ggsave("result/Fig1c.jpg", p, width = 7, height = 12, dpi = 600)
ggsave("result/Fig1c.pdf", p, width = 7, height = 12, dpi = 600)

library(ggplot2)
library(ggridges)
library(dplyr)
library(ggsignif)
library(tidyverse)

df$cluster <- str_extract(df[,1], "^[^.]+")

df_count$Database <- factor(df_count$Database.type,
                            levels = c("The Aral Sea Basin", "Public"))
df_count$Quality <- factor(df_count$Quality.level,
                           levels = c("Medium-quality", "High-quality"))

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ...,
                              draw_quantiles = NULL,
                              trim = TRUE, scale = "area",
                              na.rm = TRUE, orientation = NA,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles,
                      na.rm = na.rm, orientation = orientation, ...))
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (xmax - x),
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
                                                      if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ])
                             if (nrow(newdata) > 1) {
                               newdata[c(1, nrow(newdata)), "x"] <- newdata[1, "x"]
                             }
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)),
                                                  setdiff(names(data), c("y", "x")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob)
                             } else {
                               GeomPolygon$draw_panel(newdata, ...)
                             }
                           })
library(ggplot2)
library(ggsignif)
library(plyr)

ggplot(df_count, aes(
  x = Database,
  y = MAG_count,
  fill = Quality
)) +
  geom_split_violin(alpha = 0.7, trim = FALSE, color = NA) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 2,
               position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15,
               position = position_dodge(width = 0.2)) +
  stat_compare_means(
    aes(group = Quality),
    method = "wilcox.test",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "NS")),
    label = "p.signif") +

scale_fill_manual(
  values = c("Medium-quality" = "#FDE725", "High-quality" = "#5DC863"),
  name = "Quality Level"
) +
  labs(
    x = "",
    y = "Number of MAGs per SGB"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    panel.border = element_rect(linewidth = 0.8, fill = NA, colour = "black"),
    axis.line = element_line(colour = "black")
  )
ggsave("result/Fig1d.tiff", width = 3.5, height = 4, dpi = 300)
ggsave("result/Fig1d.jpg", width = 3.5, height = 4, dpi = 300)
ggsave("result/Fig1d.pdf", width = 4, height = 4, dpi = 300)
p_e <- ggplot(df, aes(x = Genome.size, y = completeness, color = Quality.level)) +
  geom_point(alpha = 0.3, size=1.2) +
  scale_color_manual(values = c("High-quality"="#5DC863","Medium-quality"="#FDE725")) +
  labs(x="Genome.size", y="completeness") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        legend.position = "none")
p_e
ggsave("result/Fig1e.tiff",p_e, width = 3, height = 4, dpi = 300)
ggsave("result/Fig1e.jpg",p_e, width = 3, height = 4, dpi = 300)
ggsave("result/Fig1e.pdf",p_e, width = 3, height = 4, dpi = 300)
p_f <- ggplot(df, aes(x = Genome.size, y = contamination, color = Quality.level)) +
  geom_point(alpha = 0.3, size=1.2) +
  scale_color_manual(values = c("High-quality"="#5DC863","Medium-quality"="#FDE725")) +
  labs(x="Genome.size", y="contamination") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        legend.position = "none")
p_f
ggsave("result/Fig1f.tiff",p_f, width = 3, height = 4, dpi = 300)
ggsave("result/Fig1f.jpg",p_f, width = 3, height = 4, dpi = 300)
ggsave("result/Fig1f.pdf",p_f, width = 3, height = 4, dpi = 300)
p_g <- ggplot(df, aes(x = completeness, y = contamination, color = Quality.level)) +
  geom_point(alpha = 0.3, size=1.2) +
  scale_color_manual(values = c("High-quality"="#5DC863","Medium-quality"="#FDE725")) +
  labs(x="completeness", y="contamination") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.3),
        legend.position = "none")
p_g
ggsave("result/Fig1g.tiff",p_g, width = 3, height = 4, dpi = 300)
ggsave("result/Fig1g.jpg",p_g, width = 3, height = 4, dpi = 300)
ggsave("result/Fig1g.pdf",p_g, width = 3, height = 4, dpi = 300)

library(tidyverse)

bins_breaks <- c(0, 0.0001,10,20,30,40,50,60,70,80,90,100)
bins_labels <- c("0","0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")

df_stack <- df %>%
  mutate(
    bin_range = cut(
      strain_heterogeneity,
      breaks = bins_breaks,
      labels = bins_labels,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  group_by(bin_range, Quality.level) %>%
  summarise(count = n(), .groups = "drop")

df_stack$Quality.level <- factor(df_stack$Quality.level,
                                 levels = c("Medium-quality", "High-quality"))
df_stack$bin_range <- factor(df_stack$bin_range, levels = bins_labels)

p_h <- ggplot(df_stack, aes(x = bin_range, y = count, fill = Quality.level)) +
  geom_col(position = position_stack(reverse = TRUE),
           width = 0.85, color = "white", linewidth=0.25) +
  scale_fill_manual(
    values = c("High-quality"="#007020", "Medium-quality"="#a4d1bf"),
    name = "Quality Level"
  ) +
  labs(
    x = "Strain heterogeneity (%)",
    y = "SGBs"
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(linewidth = 0.6, fill = NA),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.83, 0.86),
    legend.background = element_blank()
  )+
  guides(fill = guide_legend(reverse = TRUE))
p_h
ggsave("result/Fig1h.tiff",p_h, width = 8, height = 4, dpi = 300)
ggsave("result/Fig1h.jpg",p_h, width = 8, height = 4, dpi = 300)
ggsave("result/Fig1h.pdf",p_h, width = 8, height = 4, dpi = 300)

library(tidyverse)
df <- read.csv("data/SGB_info.csv", fileEncoding = "GBK")
phylum_order <- c(
  "Pseudomonadota",
  "Actinomycetota",
  "Bacteroidota",
  "Halobacteriota",
  "Thermoplasmatota",
  "Others"
)
phylum_color <- c(
  "Pseudomonadota" = "#d8d2c9",
  "Actinomycetota" = "#b8a399",
  "Bacteroidota"   = "#7c88a0",
  "Halobacteriota" = "#5b7494",
  "Thermoplasmatota" = "#527d64",
  "Others" = "#c6c8ca"
)

p_a_data <- df %>%
  count(Phylum..gtdb., name = "n") %>%
  rename(Phylum = Phylum..gtdb.) %>%
  mutate(Phylum = case_when(
    Phylum %in% phylum_order[1:5] ~ Phylum,
    TRUE ~ "Others"
  )) %>%
  group_by(Phylum) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(percent = paste0(round(n/sum(n)*100,1),"%")) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order))

p_a <- ggplot(p_a_data, aes(x = 2, y = n, fill = Phylum)) +
  geom_col(color = "white", linewidth=0.4) +
  geom_text(aes(label = percent), position = position_stack(vjust = 0.5), size = 3.3) +
  coord_polar("y", start = 0) +
  xlim(0, 2.6) +
  scale_fill_manual(values = phylum_color) +
  labs(title = "a") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    legend.position = "bottom"
  )
ggsave("result/sFig1a.tiff",p_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig1a.jpg",p_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig1a.pdf",p_a, width = 5, height = 4, dpi = 300)
p_b_data <- df %>%
  filter(Quality.level == "High-quality") %>%
  count(Phylum..gtdb., name = "n") %>%
  rename(Phylum = Phylum..gtdb.) %>%
  mutate(Phylum = case_when(
    Phylum %in% phylum_order[1:5] ~ Phylum,
    TRUE ~ "Others"
  )) %>%
  group_by(Phylum) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(percent = paste0(round(n/sum(n)*100,1),"%")) %>%
  mutate(Phylum = factor(Phylum, levels = phylum_order))

p_b <- ggplot(p_b_data, aes(x = 2, y = n, fill = Phylum)) +
  geom_col(color = "white", linewidth=0.4) +
  geom_text(aes(label = percent), position = position_stack(vjust = 0.5), size = 3.3) +
  coord_polar("y", start = 0) +
  xlim(0, 2.6) +
  scale_fill_manual(values = phylum_color) +
  labs(title = "b") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    legend.position = "bottom"
  )
ggsave("result/sFig1b.tiff",p_b, width = 5, height = 4, dpi = 300)
ggsave("result/sFig1b.jpg",p_b, width = 5, height = 4, dpi = 300)
ggsave("result/sFig1b.pdf",p_b, width = 5, height = 4, dpi = 300)

library(tidyverse)
library(FSA)
df <- read.csv("data/SGB_info.csv", header = TRUE)
df <- df %>%
  mutate(
    Kingdom = str_trim(Kingdom..gtdb.),
    Quality.level = str_trim(Quality.level)
  ) %>%
  mutate(
    group_type = case_when(
      Kingdom..gtdb. == "Archaea" & Quality.level == "High-quality" ~ "Archaea.High",
      Kingdom..gtdb. == "Archaea" & Quality.level == "Medium-quality" ~ "Archaea.Medium",
      Kingdom..gtdb. == "Bacteria" & Quality.level == "High-quality" ~ "Bacteria.High",
      Kingdom..gtdb. == "Bacteria" & Quality.level == "Medium-quality" ~ "Bacteria.Medium"
    )
  ) %>% filter(!is.na(group_type))

group_order <- c("Archaea.High","Bacteria.High","Archaea.Medium","Bacteria.Medium")
df$group_type <- factor(df$group_type, levels = group_order)

get_signif_label <- function(data, yvar){
  formula_str <- as.formula(paste0(yvar, " ~ group_type"))
  dt <- dunnTest(formula_str, data=data, method="fdr")
  res <- dt$res

  cat("\n=====", yvar, " FDR-adjusted P values =====\n")
  print(res[,c("Comparison","P.adj")])

  res$signif <- case_when(
    res$P.adj < 0.0001 ~ "****",
    res$P.adj < 0.001 ~ "***",
    res$P.adj < 0.01 ~ "**",
    res$P.adj < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  return(res)
}

get_signif_label <- function(data, yvar){
  formula_str <- as.formula(paste0(yvar, " ~ group_type"))
  dt <- dunnTest(formula_str, data=data, method="bh") # 重点修改这里！
  res <- dt$res

  cat("\n=====", yvar, " FDR(BH)-adjusted P values =====\n")
  print(res[,c("Comparison","P.adj")])

  res$signif <- case_when(
    res$P.adj < 0.0001 ~ "****",
    res$P.adj < 0.001 ~ "***",
    res$P.adj < 0.01 ~ "**",
    res$P.adj < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  return(res)
}

sig_gc <- get_signif_label(df, "GC")
sig_gs <- get_signif_label(df, "Genome.size")

df_long <- df %>%
  pivot_longer(
    cols = c("GC", "Genome.size"),
    names_to = "variable",
    values_to = "Value"
  ) %>%
  mutate(variable = factor(variable,
                           levels = c("GC","Genome.size"),
                           labels = c("GC content (%)","Genome size")))

fill_col <- c("High-quality"="#5DC863","Medium-quality"="#FDE725")
p_c <- ggplot(df_long, aes(x = group_type, y = Value)) +
  geom_boxplot(aes(fill=Quality.level), width=0.4, outlier.shape=NA) +
  geom_jitter(width = 0.12, size=1, alpha=0.2, color="black") +
  scale_fill_manual(values = fill_col) +
  facet_wrap(.~variable, scales="free_y", nrow=1) +
  labs(x="", y="Value") +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle=35, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_c
ggsave("result/sFig1c.tiff",p_c, width = 8, height = 8, dpi = 300)
ggsave("result/sFig1c.jpg",p_c, width = 8, height = 8, dpi = 300)
ggsave("result/sFig1c.pdf",p_c, width = 8, height = 8, dpi = 300)

library(tidyverse)
library(scales)
library(patchwork)
phylum_colors <- c(
  "Others" = "#F2F2E9", "Bacteroidota" = "#A3C5A8",
  "Gemmatimonadota" = "#DCD6E8", "Pseudomonadota" = "#5AB4AC",
  "Acidobacteriota" = "#2D728F", "Chloroflexota" = "#E8D8D3",
  "Halobacteriota" = "#9AB87A", "Thermoplasmatota" = "#F29468",
  "Actinomycetota" = "#D85A5A", "Cyanobacteriota" = "#D89A78",
  "Myxococcota" = "#995D8D", "Thermoproteota" = "#A3D4E8",
  "Bacillota" = "#3D7C47", "Desulfobacterota" = "#F2D06B",
  "Planctomycetota" = "#2D558F", "Verrucomicrobiota" = "#C8C8C8"
)

phylum_all <- names(phylum_colors)
donut_raw <- df %>%
  group_by(Phylum..gtdb.,Novel_GTDB) %>%
  rename(Phylum = Phylum..gtdb.) %>%
  mutate(Phylum = ifelse(Phylum %in% phylum_all, Phylum, "Others")) %>%
  summarise(Count = n(), .groups = "drop")

full_grid <- expand.grid(
  Phylum.gtdb. = phylum_all,
  Novel_GTDB = c("Novel", "Defined")
)
colnames(donut_raw)
donut_data <- donut_raw %>%
  complete(
    Phylum = phylum_all,
    Novel_GTDB = c("Novel","Defined"),
    fill = list(Count = 0)
  )

donut_data <- donut_data %>%
  group_by(Novel_GTDB) %>%
  mutate(Percent = Count/sum(Count)*100) %>%
  ungroup()

p_2a <- ggplot(donut_data, aes(x = 2, y = Percent, fill = Phylum)) +
  geom_col(width = 1, color = NA) +
  coord_polar("y", start = 0) +
  facet_wrap(~Novel_GTDB, nrow = 1) +
  xlim(1, 2.5) +
  scale_fill_manual(values = phylum_colors) +
  theme_void() +
  theme(
    strip.text = element_text(size = 16),
    legend.position = "bottom"
  )

p_2a
ggsave("result/sFig2a.tiff",p_2a, width = 6, height = 4, dpi = 300)
ggsave("result/sFig2a.jpg",p_2a, width = 6, height = 4, dpi = 300)
ggsave("result/sFig2a.pdf",p_2a, width = 6, height = 4, dpi = 300)

class_colors <- c(
  "Acidimicrobiia" = "#E8E9CC",
  "Actinomycetes" = "#247B8C",
  "Alphaproteobacteria" = "#D04A42",
  "Anaerolineae" = "#247037",
  "Bacilli" = "#88B499",
  "Bacteroidia" = "#E8D9D4",
  "Gammaproteobacteria" = "#E09566",
  "Gemmatimonadetes" = "#F2D058",
  "Halobacteria" = "#DCD6EC",
  "Nitriliruptoria" = "#88B040",
  "Nitrososphaeria" = "#80398C",
  "Planctomycetia" = "#2D62AF",
  "Polyangia" = "#48AA90",
  "Rhodothermia" = "#F08950",
  "Thermoanaerobaculia" = "#A6C8E0",
  "Others" = "#F2F2E9"
)

class_all <- names(class_colors)
donut_raw <- df %>%
  rename(Class = Class..gtdb.) %>%
  mutate(Class = ifelse(Class %in% class_all, Class, "Others")) %>%
  group_by(Class, Novel_GTDB) %>%
  summarise(Count = n(), .groups = "drop")

donut_data <- donut_raw %>%
  complete(
    Class = class_all,
    Novel_GTDB = c("Novel","Defined"),
    fill = list(Count = 0)
  )

donut_data <- donut_data %>%
  group_by(Novel_GTDB) %>%
  mutate(Percent = Count/sum(Count)*100) %>%
  ungroup()

p_2b <- ggplot(donut_data, aes(x = 2, y = Percent, fill = Class)) +
  geom_col(width = 1, color = NA) +
  coord_polar("y", start = 0) +
  facet_wrap(~Novel_GTDB, nrow = 1) +
  xlim(1, 2.5) +
  scale_fill_manual(values = class_colors, name = "Class") +
  guides(fill = guide_legend(nrow = 4)) +
  theme_void() +
  theme(
    strip.text = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

p_2b
ggsave("result/sFig2b.tiff",p_2b, width = 6, height = 4, dpi = 300)
ggsave("result/sFig2b.jpg",p_2b, width = 6, height = 4, dpi = 300)
ggsave("result/sFig2b.pdf",p_2b, width = 6, height = 4, dpi = 300)
library(tidyverse)
library(patchwork)
library(readr)
library(readxl)
tax <- read.csv("data/SGB_taxonomy.csv", fileEncoding = "GBK")
info <- read.csv("data/SGB_info.csv", fileEncoding = "GBK")
dat <- inner_join(tax, info)
colnames(tax)
colnames(info)
dat <- dat %>%
  mutate(
    RED = `dRep.Score` / 100,

    Novelty = case_when(
      `Species..gtdb.` != "Unclassified" ~ "SGB",
      `Genus..gtdb.` != "Unclassified" ~ "uGGB",
      `Family..gtdb.` != "Unclassified" ~ "uFGB",
      `Order..gtdb.` != "Unclassified" ~ "uOGB",
      TRUE ~ "uSGB"
    ))
keep_phylum <- dat %>%
  count(`Phylum..gtdb.`) %>%
  filter(n >= 3) %>%
  pull(`Phylum..gtdb.`)

dat_c <- dat %>% filter(`Phylum..gtdb.` %in% keep_phylum)
p_c <- ggplot(dat_c, aes(x = `Phylum..gtdb.`, y = `size`)) +
  geom_boxplot(fill="#87CEEB", outlier.shape = NA)+
  geom_jitter(width=0.1, size=0.5)+
  labs(y = "Genome Size (bp)", x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank())
p_c
ggsave("result/sFig2c.tiff",p_c, width = 12, height = 4, dpi = 300)
ggsave("result/sFig2c.jpg",p_c, width = 12, height = 4, dpi = 300)
ggsave("result/sFig2c.pdf",p_c, width = 12, height = 4, dpi = 300)
target_phylum = c("Halobacteriota", "Nanohaloarchaeota", "Thermoplasmatota", "Thermoproteota","Nanoarchaeota")
dat_d <- dat %>%
  filter(`Phylum..gtdb.` %in% target_phylum) %>%
  mutate(`Phylum..gtdb.` = factor(`Phylum..gtdb.`, levels = target_phylum))

p_d <- ggplot(dat_d, aes(y=`Phylum..gtdb.`, x = `size`))+
  geom_boxplot(fill="#87CEEB", outlier.shape=NA)+
  geom_jitter(height=0.2, size=0.8)+
  labs(x = "Genome Size (bp)", y="")+
  theme_bw()+
  theme(axis.text.y = element_text(angle = 60, hjust = 1),panel.grid = element_blank())
p_d
ggsave("result/sFig2d.tiff",p_d, width = 3, height = 6, dpi = 300)
ggsave("result/sFig2d.jpg",p_d, width = 3, height = 6, dpi = 300)
ggsave("result/sFig2d.pdf",p_d, width = 3, height = 6, dpi = 300)

library(tidyverse)

# 读取数据
df1 <- read.csv("F:/Saline-alkaline-Microbiome/Metagenome/MAG_stat/data/SGB_info.csv", row.names = NULL, fileEncoding = "GBK")
df_processed <- df1 %>%
  mutate(
    red_value = ifelse(red_value == "N/A", NA, red_value),
    red_value = as.numeric(red_value),
    group = case_when(
      # 从门到种，由高到低判断
      Phylum..gtdb. == "Unclassified" ~ "uPGB",
      Class..gtdb. == "Unclassified" ~ "uCGB",
      Order..gtdb. == "Unclassified" ~ "uOGB",
      Family..gtdb. == "Unclassified" ~ "uFGB",
      Genus..gtdb. == "Unclassified" ~ "uGGB",
      Species..gtdb. == "Unclassified" ~ "uSGB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group), !is.na(red_value))

# 查看4组数量
table(df_processed$group)

# 绘图：4个面板
p_e <- ggplot(df_processed, aes(x = red_value)) +
  geom_histogram(aes(fill = group), bins = 50, color = NA) +
  facet_wrap(~group, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = c(
    "uFGB" = "#2ca02c",
    "uGGB" = "#1f77b4",
    "uOGB" = "#9467bd",
    "uSGB" = "#d62728"
  )) +
  labs(x = "RED value", y = "Count") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray80"),
    panel.grid = element_blank()
  )
print(p_e)
ggsave("result/sFig2e.tiff",p_e, width = 8, height = 6, dpi = 300)
ggsave("result/sFig2e.jpg",p_e, width = 8, height = 6, dpi = 300)
ggsave("result/sFig2e.pdf",p_e, width = 8, height = 6, dpi = 300)
library(patchwork)

final_fig <- (p_2a + p_2b) / p_c / (p_d + p_e) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(plot.tag.position = c(0.04,0.96))

final_fig
ggsave("result/sFig2.tiff",final_fig, width = 14, height = 20, dpi = 300)
ggsave("result/sFig2.jpg",final_fig, width = 14, height = 20, dpi = 300)
ggsave("result/sFig2.pdf",final_fig, width = 14, height = 20, dpi = 300)
library(tidyverse)
library(patchwork)

#==================== 读取数据 ====================
df_raw <- read.csv("data/SGB_info.csv", header = TRUE, fileEncoding = "GBK") # 修改成你的文件路径

library(tidyverse)
library(patchwork)

# ===================== 定义目标分类列表 =====================
# 门水平 保留15个phylum
target_phylum <- c(
  "Acidobacteriota",
  "Actinomycetota",
  "Bacillota",
  "Bacteroidota",
  "Chloroflexota",
  "Cyanobacteriota",
  "Desulfobacterota",
  "Gemmatimonadota",
  "Halobacteriota",
  "Myxococcota",
  "Planctomycetota",
  "Pseudomonadota",
  "Thermoplasmatota",
  "Thermoproteota",
  "Verrucomicrobiota"
)

# 纲水平，图b里面展示的所有Class
target_class <- c(
  "Acidimicrobiia",
  "Actinomycetes",
  "Alphaproteobacteria",
  "Anaerolineae",
  "Bacilli",
  "Bacteroidia",
  "Gammaproteobacteria",
  "Gemmatimonadetes",
  "Halobacteria",
  "Nitriliruptoria",
  "Nitrososphaeria",
  "Planctomycetia",
  "Polyangia",
  "Rhodothermia",
  "Thermoanaerobaculia"
)

# ===================== 合并分类函数：不在名单→Others =====================
merge_taxon <- function(tbl, tax_col, target_list){
  tbl %>%
    count(!!sym(tax_col), name="count") %>%
    rename(raw_taxon = 1) %>%
    mutate(taxon = ifelse(raw_taxon %in% target_list, raw_taxon, "Others")) %>%
    group_by(taxon) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    # 设置因子，固定顺序
    mutate(taxon = factor(taxon, levels = c(target_list, "Others"))) %>%
    mutate(is_other = taxon == "Others") %>%
    arrange(is_other, desc(count)) %>%
    select(-is_other)
}

# 拆分Presence / Unclassified（按Species是否Unclassified）
df_pres <- df_raw %>% filter(Species..gtdb. != "Unclassified")
df_uncl <- df_raw %>% filter(Species..gtdb. == "Unclassified")
df_raw
df_aral <- df_raw %>% filter(Database.type == "The Aral Sea Basin")
df_aral
df_pub <- df_raw %>% filter(Database.type == "Public")

# 门水平数据
ph_pres <- merge_taxon(df_pres, "Phylum..gtdb.", target_phylum)
ph_uncl <- merge_taxon(df_uncl, "Phylum..gtdb.", target_phylum)

ph_aral <- merge_taxon(df_aral, "Phylum..gtdb.", target_phylum)
ph_pub <- merge_taxon(df_pub, "Phylum..gtdb.", target_phylum)

# 纲水平数据
cl_pres <- merge_taxon(df_pres, "Class..gtdb.", target_class)
cl_uncl <- merge_taxon(df_uncl, "Class..gtdb.", target_class)

cl_aral <- merge_taxon(df_aral, "Class..gtdb.", target_class)
cl_pub <- merge_taxon(df_pub, "Class..gtdb.", target_class)

pal_phylum <- c(
  "Pseudomonadota" = "#5AB4AC",
  "Actinomycetota" = "#D85A5A",
  "Bacteroidota" = "#A3C5A8",
  "Acidobacteriota" = "#2D728F",
  "Chloroflexota" = "#E8D8D3",
  "Gemmatimonadota" = "#DCD6E8",
  "Halobacteriota" = "#9AB87A",
  "Planctomycetota" = "#2D558F",
  "Desulfobacterota" = "#F2D06B",
  "Myxococcota" = "#995D8D",
  "Thermoplasmatota" = "#F29468",
  "Cyanobacteriota" = "#D89A78",
  "Bacillota" = "#3D7C47",
  "Thermoproteota" = "#A3D4E8",
  "Verrucomicrobiota" = "#C8C8C8",
  "Others" = "#F2F2E9"
)

pal_class <- c(
  "Acidimicrobiia"="#A3C5A8",
  "Actinomycetes"="#DCD6E8",
  "Alphaproteobacteria"="#5AB4AC",
  "Anaerolineae"="#2D728F",
  "Bacilli"="#E8D8D3",
  "Bacteroidia"="#9AB87A",
  "Gammaproteobacteria"="#F29468",
  "Gemmatimonadetes"="#D85A5A",
  "Halobacteria"="#D89A78",
  "Nitriliruptoria"="#995D8D",
  "Nitrososphaeria"="#A3D4E8",
  "Planctomycetia"="#3D7C47",
  "Polyangia"="#F2D06B",
  "Rhodothermia"="#2D558F",
  "Thermoanaerobaculia"="#C8C8C8",
  "Others"="#F6F6F0"
)

# ===================== 绘图函数：甜甜圈图，子图无图例 =====================
donut_plot <- function(dat, title, pal){
  ggplot(dat, aes(x=2, y=count, fill=taxon)) +
    geom_col(color="white", linewidth=0.2) +
    xlim(1, 2.8) +
    coord_polar("y", start = pi/2) +
    scale_fill_manual(values = pal, drop=TRUE) +
    labs(title=title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=16),
      legend.position = "none" #子图内部不要图例
    )
}

# 绘制4张子图
p_ph_aral <- donut_plot(ph_aral, "Aral", pal_phylum)
p_ph_pub <- donut_plot(ph_pub, "Public", pal_phylum)
p_ph_aral
p_ph_pub
p_cl_aral <- donut_plot(cl_aral, "Aral", pal_class)
p_cl_pub <- donut_plot(cl_pub, "Public", pal_class)

# ===================== 提取图例（单独提取，放在底部） =====================
get_legend_plot <- function(dat, pal, ncol=4){
  ggplot(dat, aes(x=1, y=count, fill=taxon)) +
    geom_col() +
    scale_fill_manual(values=pal, drop=TRUE) +
    guides(fill = guide_legend(ncol = ncol)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

leg_ph <- get_legend_plot(ph_aral, pal_phylum, ncol=4)
leg_cl <- get_legend_plot(cl_aral, pal_class, ncol=4)
leg_ph
leg_cl
# ===================== 拼图：a图(门)一行2图 + 底部图例；b图(纲)一行2图 + 底部图例 =====================
plot_a <- (p_ph_aral + p_ph_pub) / leg_ph + plot_layout(heights = c(10, 1))
plot_b <- (p_cl_aral + p_cl_pub) / leg_cl + plot_layout(heights = c(10, 1))
plot_a
plot_b

library(ggpubr)
final_fig <- ggarrange(
  plot_a, plot_b,
  ncol = 1, nrow = 2,
  labels = c("a","b"),
  label.x = 0.04, label.y = 0.96,
  common.legend = FALSE
)
final_fig

ggsave("result/sFig3a.tiff",plot_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3a.jpg",plot_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3a.pdf",plot_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3b.tiff",plot_b, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3b.jpg",plot_b, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3b.pdf",plot_b, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3c.tiff",plot_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3c.jpg",plot_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3c.pdf",plot_a, width = 5, height = 4, dpi = 300)
ggsave("result/sFig3.tiff",final_fig, width = 6, height = 10, dpi = 300)
ggsave("result/sFig3.jpg",final_fig, width = 6, height = 10, dpi = 300)
ggsave("result/sFig3.pdf",final_fig, width = 6, height = 10, dpi = 300)

# 安装包（只跑一次）
#BiocManager::install(c("ggtree","ggtreeExtra","treeio"))
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(dplyr)
library(cowplot)
library(grid)
library(ggplot2)
library(ggnewscale)

#==================== 1 读入数据 ====================
tree <- read.tree("data/gtdbtk_iqtree2_result.treefile")   # SGB进化树newick
metadata <- read.csv("data/SGB_info.csv", row.names = NULL, fileEncoding = "GBK")
# metadata列：GenomeID, Group, Taxonomy, Species_size, Novel_GTDB

novel_color <- c("Novel"="#9acd00",p_e <- ggtree(tree,layout="circular",size=0.2) %<+% metadata +
  # 内层：门颜色条 tile
  geom_tile(aes(x=x+0.3,y=y,fill=phylum), inherit.aes=F, width=0.4)+
  # ========== 外圈防御系统条形，geom_fruit ==========
  geom_fruit(
    data=tree_def_anno,
    geom=geom_col,
    mapping=aes(x=defense_total), # 绑定树之后，不需要写y！！
    fill="#ffbc69",
    offset=1.5,
    pwidth=0.3
  )+
  scale_fill_manual(values = c(phy_color, "Defined"="#ffffff"))+
  labs(fill="Phylum")+
  theme(legend.position="bottom"))
group_color <- c(
  "Published genomes"="#c5e0e8",
  "Unclassified public database SGBs"="#d9d2e9",
  "Unclassified Aral Sea Basin SGBs"="#b4a7d6"
)

# ---------------------- 1. 定义颜色 + 预处理metadata（必须先跑！） ----------------------
tax_color <- c(
  "Acidobacteriota"="#b6d7a8",
  "Actinomycetota"="#b4a79b",
  "Bacillota"="#f4cccc",
  "Bacteroidota"="#6fa8dc",
  "Chloroflexota"="#e5efd9",
  "Cyanobacteriota"="#c27ba0",
  "Desulfobacterota"="#f8cbad",
  "Gemmatimonadota"="#76a5af",
  "Myxococcota"="#e6b89c",
  "Others"="#dddddd",
  "Planctomycetota"="#ea9999",
  "Pseudomonadota"="#d0c9c0",
  "Verrucomicrobiota"="#996644"
)

keep_phyla <- setdiff(names(tax_color), "Others")
metadata <- metadata %>%
  mutate(
    `Phylum..gtdb.` = ifelse(`Phylum..gtdb.` %in% keep_phyla, `Phylum..gtdb.`, "Others")
  )

#校验：确认是否生成Others
unique(metadata$`Phylum..gtdb.`)

# ---------------------- 2. 你原来的树和第一圈 geom_aline（Group） ----------------------
tree.plot <- ggtree(tree, layout = "fan", right = TRUE, open.angle = 15, size = 0.1) %<+% metadata +
  geom_aline(aes(color = Group), linetype = "solid", size = 0.2, show.legend = TRUE) +
  scale_color_manual(values = group_color, name = "Group") +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
    )
  ) +
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.5, "mm")
  )
tree.plot
# ---------------------- 3. 在后面叠加第二圈 Phylum tile【关键！new_scale_fill()】 ----------------------
main.plot <- tree.plot +
  new_scale_fill() +  # 新增fill比例尺，解决多层叠加冲突！
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = `Phylum..gtdb.`), # y这里用label！！不是User_genome！！
    width = 0.06,
    pwidth = 0.05
  ) +
  scale_fill_manual(values = tax_color, name = "Taxonomy") +
  guides(fill = guide_legend(
    title.position = "top",
    title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
  ))

main.plot
colnames(metadata)
main.plot.temp <- main.plot +
  new_scale_fill() +
  geom_fruit(
    geom = geom_bar,
    stat = "identity",
    width = 0.7,
    aes(y = label, x = log10(cluster_members)),
    pwidth = 0.25
  ) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    width = 0.3,
    aes(y = label, fill = Novel_GTDB),
    pwidth = 0.05
  ) +
  scale_fill_manual(values = novel_color, name = "Novelty")
main.plot.temp
ggsave("result/Fig2a.tiff",main.plot.temp, width = 15, height = 10, dpi = 300)
ggsave("result/Fig2a.jpg",main.plot.temp, width = 15, height = 10, dpi = 300)
ggsave("result/Fig2a.pdf",main.plot.temp, width = 15, height = 10, dpi = 300)

tree <- read.tree("data/gtdbtk.unrooted2.tree")   # SGB进化树newick
metadata <- read.csv("data/SGB_info.csv", row.names = NULL)
# metadata列：GenomeID, Group, Taxonomy, Species_size, Novel_GTDB

novel_color <- c("Novel"="#9acd00",
                 "Defined"="#ffffff")

# ---------------------- 1. 定义颜色 + 预处理metadata（必须先跑！） ----------------------
tax_color <- c(
  "Halobacteriota" = "#4a5f76",
  "Thermoplasmatota" = "#4c6b57",
  "Thermoproteota" = "#d9b38c",
  "Nanoarchaeota"="#b0cfc7",
  "Nanohaloarchaeota"="#aaabc2",
  "SpSt-1190"="#cae7df"
)

keep_phyla <- setdiff(names(tax_color), "Others")
metadata <- metadata %>%
  mutate(
    `Phylum..gtdb.` = ifelse(`Phylum..gtdb.` %in% keep_phyla, `Phylum..gtdb.`, "Others")
  )

#校验：确认是否生成Others
unique(metadata$`Phylum..gtdb.`)

# ---------------------- 2. 你原来的树和第一圈 geom_aline（Group） ----------------------
tree.plot <- ggtree(tree, layout = "fan", right = TRUE, open.angle = 15, size = 0.1) %<+% metadata +
  geom_aline(aes(color = Group), linetype = "solid", size = 0.2, show.legend = TRUE) +
  scale_color_manual(values = group_color, name = "Group") +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
    )
  ) +
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.5, "mm")
  )
tree.plot
# ---------------------- 3. 在后面叠加第二圈 Phylum tile【关键！new_scale_fill()】 ----------------------
main.plot <- tree.plot +
  new_scale_fill() +  # 新增fill比例尺，解决多层叠加冲突！
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = `Phylum..gtdb.`), # y这里用label！！不是User_genome！！
    width = 0.06,
    pwidth = 0.05
  ) +
  scale_fill_manual(values = tax_color, name = "Taxonomy") +
  guides(fill = guide_legend(
    title.position = "top",
    title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
  ))

main.plot
colnames(metadata)
main.plot.temp <- main.plot +
  new_scale_fill() +
  geom_fruit(
    geom = geom_bar,
    stat = "identity",
    width = 0.7,
    aes(y = label, x = log10(cluster_members)),
    pwidth = 0.25
  ) +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    width = 0.1,
    aes(y = label, fill = Novel_GTDB),
    pwidth = 0.05
  ) +
  scale_fill_manual(values = novel_color, name = "Novelty")
main.plot.temp
ggsave("result/Fig2c.tiff",main.plot.temp, width = 10, height = 10, dpi = 300)
ggsave("result/Fig2c.jpg",main.plot.temp, width = 10, height = 10, dpi = 300)
ggsave("result/Fig2c.pdf",main.plot.temp, width = 10, height = 10, dpi = 300)

library(tidyverse)
library(patchwork)

df_bac_raw <- read.csv("data/Bacteria.PD.csv") %>%
  mutate(
    across(
      c(The_Aral_Sea_novelty, Public_novelty, Published),
      ~suppressWarnings(as.numeric(.)) %>% replace_na(0)
    )
  )

df_bac <- df_bac_raw %>%
  pivot_longer(
    cols = c(The_Aral_Sea_novelty, Public_novelty, Published),
    names_to = "Category",
    values_to = "PD"
  ) %>%
  # 关键：设置因子顺序，两个紫色在前，保证堆叠时相邻
  mutate(Category = factor(Category, levels = c("The_Aral_Sea_novelty", "Public_novelty", "Published")))

# 计算百分比
total_bac_pd <- sum(df_bac$PD)
df_bac <- df_bac %>% mutate(PD_percent = PD / total_bac_pd * 100)

# 按总PD排序门类
phylum_total <- df_bac %>%
  group_by(Bacteria_Phylum) %>%
  summarise(sumPD = sum(PD), .groups = "drop") %>%
  arrange(sumPD)
df_bac$Bacteria_Phylum <- factor(df_bac$Bacteria_Phylum, levels = phylum_total$Bacteria_Phylum)

# 配色（完全按你的要求）
color_map <- c(
  "The_Aral_Sea_novelty" = "#8c96c6",  # 深紫
  "Public_novelty" = "#b4b9dd",        # 浅紫
  "Published" = "#b6e2cc"              # 浅绿
)

# 绘制细菌主图
plot_b <- ggplot(df_bac, aes(y = Bacteria_Phylum, x = PD_percent, fill = Category)) +
  geom_col(position = position_stack(), width = 0.7) +
  scale_fill_manual(values = color_map) +
  labs(x = "Phylogenetic contribution of each phylum (%)", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 9)
  )

plot_b
# ====== 2. 古菌部分 ======
df_arch_raw <- read.csv("data/Archaea.PD.csv") %>%
  mutate(
    across(
      c(The_Aral_Sea_novelty, Public_novelty, Published),
      ~suppressWarnings(as.numeric(.)) %>% replace_na(0)
    )
  )

df_arch <- df_arch_raw %>%
  pivot_longer(
    cols = c(The_Aral_Sea_novelty, Public_novelty, Published),
    names_to = "Category",
    values_to = "PD"
  ) %>%
  # 同样设置因子顺序，保证和细菌图一致的堆叠逻辑
  mutate(Category = factor(Category, levels = c("The_Aral_Sea_novelty", "Public_novelty", "Published")))

total_arch_pd <- sum(df_arch$PD)
df_arch <- df_arch %>% mutate(PD_percent = PD / total_arch_pd * 100)

arch_total <- df_arch %>%
  group_by(Archaea_phylum) %>%
  summarise(sumPD = sum(PD), .groups = "drop") %>%
  arrange(sumPD)
df_arch$Archaea_phylum <- factor(df_arch$Archaea_phylum, levels = arch_total$Archaea_phylum)

# 绘制古菌内嵌小图
plot_arch <- ggplot(df_arch, aes(y = Archaea_phylum, x = PD_percent, fill = Category)) +
  geom_col(position = position_stack(), width = 0.7) +
  scale_fill_manual(values = color_map) +
  labs(x = "Phylogenetic contribution (%)", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 9)
  )
plot_arch
# ===================== 3. 拼图（主图+内嵌小图） =====================
final_plot <- plot_b + inset_element(plot_arch, left = 0.3, bottom = 0.1, right = 0.98, top = 0.62)

# 显示最终图
final_plot

ggsave("result/Fig2b.tiff",final_plot, width = 12, height = 10, dpi = 300)
ggsave("result/Fig2b.jpg",final_plot, width = 12, height = 10, dpi = 300)
ggsave("result/Fig2b.pdf",final_plot, width = 12, height = 10, dpi = 300)

#安装包（第一次运行执行）
#BiocManager::install("ComplexHeatmap")
#install.packages("circlize")
library(ComplexHeatmap)
library(circlize)

library(tidyverse)
library(readxl)

# 读取excel
df <- read.csv("data/ko_abundance_bacteria.csv", row.names = 1)

# 1. GeneRatio字符串转数值（"109/1290" → 0.0845）
df <- df %>%
  separate(GeneRatio, into = c("num", "den"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = num / den,
         neg_log_padj = -log10(p.adjust))

# 绘图
p<- ggplot(df, aes(x = GeneRatio_val, y = Description)) +
  geom_point(aes(size = Count, fill = neg_log_padj),
             shape = 21, stroke = 0.3) +
  scale_fill_gradient(low = "#f7fcb9", high = "#de2d26") +
  scale_size(range = c(2, 10)) +
  scale_y_discrete(expand = expansion(add = c(1.2,1.2))) +
  scale_x_continuous(expand = expansion(add = c(0.01,0.01))) +
  labs(
    x = "GeneRatio",
    y = "KEGG Pathway",
    fill = "-log10(p.adjust)",
    size = "Count"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t=10, r=30, b=50, l=10, unit="pt"),
    panel.grid = element_blank()
  )
p
ggsave("result/sFig4a.tiff",p, width = 10, height = 10, dpi = 300)
ggsave("result/sFig4a.jpg",p, width = 10, height = 10, dpi = 300)
ggsave("result/sFig4a.pdf",p, width = 10, height = 10, dpi = 300)

df <- read.csv("data/ko_abundance_archaea.csv", row.names = 1)

# 1. GeneRatio字符串转数值（"109/1290" → 0.0845）
df <- df %>%
  separate(GeneRatio, into = c("num", "den"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_val = num / den,
         neg_log_padj = -log10(p.adjust))

# 绘图
p1<- ggplot(df, aes(x = GeneRatio_val, y = Description)) +
  geom_point(aes(size = Count, fill = neg_log_padj),
             shape = 21, stroke = 0.3) +
  scale_fill_gradient(low = "#f7fcb9", high = "#de2d26") +
  scale_size(range = c(2, 10)) +
  scale_y_discrete(expand = expansion(add = c(1.2,1.2))) + 
  scale_x_continuous(expand = expansion(add = c(0.015,0.015))) +
  labs(
    x = "GeneRatio",
    y = "KEGG Pathway",
    fill = "-log10(p.adjust)",
    size = "Count"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(t=10, r=30, b=50, l=10, unit="pt"),
    panel.grid = element_blank()
  )
p1
ggsave("result/sFig4b.tiff",p1, width = 9, height = 10, dpi = 300)
ggsave("result/sFig4b.jpg",p1, width = 9, height = 10, dpi = 300)
ggsave("result/sFig4b.pdf",p1, width = 9, height = 10, dpi = 300)
p_combine <- p + p1 + plot_layout(widths = c(6,4)) + plot_annotation(tag_levels = "a")
p_combine
ggsave("result/sFig4.tiff",p_combine, width = 15, height = 10, dpi = 300)
ggsave("result/sFig4.jpg",p_combine, width = 15, height = 10, dpi = 300)
ggsave("result/sFig4.pdf",p_combine, width = 15, height = 10, dpi = 300)


library(tidyverse)

# ====================== 1. 读取数据+基础配置 ======================
anno_df <- read.csv("data/Bacteria.Eggnog.csv", stringsAsFactors = F)
meta <- read.csv("data/SGB_taxonomy.csv", stringsAsFactors = F)
# ID列就是基因组/bin名称！
genome_col <- "ID"  

# 8个目标功能定义
assign_function <- function(desc){
  if(is.na(desc)) return(NA)
  desc_low <- tolower(desc)
  if(str_detect(desc_low, "phosphorus|phosphate|phosphonate")) return("Phosphorus")
  if(str_detect(desc_low, "siderophore")) return("Siderophore")
  if(str_detect(desc_low, "indole")) return("Indole")
  if(str_detect(desc_low, "cytokinin")) return("Cytokinin")
  if(str_detect(desc_low, "glutathione")) return("Glutathione")
  if(str_detect(desc_low, "trehalose")) return("Trehalose")
  if(str_detect(desc_low, "glycine betaine|glycine-betaine|betaine")) return("Glycine betaine")
  if(str_detect(desc_low, "pqq")) return("PQQ")
  return(NA)
}

# ====================== 2. 清洗+整理功能-KO映射表 ======================
anno_clean <- anno_df %>%
  mutate(Function = map_chr(Description, assign_function)) %>%
  filter(!is.na(Function), KEGG_ko != "-", !is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_remove(KEGG_ko, "^ko:")) %>%
  distinct(Function, KEGG_ko, .data[[genome_col]])
func_ko_set <- anno_clean %>%
  distinct(Function, KEGG_ko) %>%
  group_by(Function) %>%
  summarise(all_ko = list(KEGG_ko), ko_num = n()) %>%
  ungroup()

func_ko_set
# ====================== 替换这段：逐基因组判定PGP功能存在性 ======================
min_count <- 1  # 至少检出1个该功能KO就算存在

genome_func_matrix <- anno_clean %>%
  group_by(ID) %>%
  summarise(genome_ko = list(unique(KEGG_ko))) %>%
  rowwise() %>%
  mutate(
    Phosphorus = {
      ko_target <- func_ko_set %>% filter(Function == "Phosphorus") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Siderophore = {
      ko_target <- func_ko_set %>% filter(Function == "Siderophore") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Indole = {
      ko_target <- func_ko_set %>% filter(Function == "Indole") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Cytokinin = {
      ko_target <- func_ko_set %>% filter(Function == "Cytokinin") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Glutathione = {
      ko_target <- func_ko_set %>% filter(Function == "Glutathione") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Trehalose = {
      ko_target <- func_ko_set %>% filter(Function == "Trehalose") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Glycine_betaine = {
      ko_target <- func_ko_set %>% filter(Function == "Glycine betaine") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    PQQ = {
      ko_target <- func_ko_set %>% filter(Function == "PQQ") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    }
  ) %>%
  ungroup()

# 保存0/1矩阵（去掉list列genome_ko！！）
genome_func_export <- genome_func_matrix %>% select(-genome_ko)

write.csv(genome_func_export, "result/PGP_function.bacteria_01_matrix.csv", row.names = F)

# ====================== 统计检出比例+绘图（直接跑） ======================
total_genome <- nrow(genome_func_export)

func_summary <- genome_func_export %>%
  summarise(across(!ID, sum)) %>%  # !ID 表示除了ID以外所有列
  pivot_longer(everything(), names_to = "Function", values_to = "Present_genome") %>%
  mutate(
    Ratio = Present_genome / total_genome * 100,
    Function = str_replace(Function, "_", " ")
  ) %>%
  arrange(desc(Ratio))

print(func_summary)
write.csv(func_summary, "result/PGP_function_bacteria_ratio.csv", row.names = F)

# 绘图
ggplot(func_summary, aes(x = reorder(Function, Ratio), y = Ratio, fill = Function)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(Ratio,1), "%")), hjust = -0.2, size = 4) +
  coord_flip() +
  labs(
    x = "PGP Functional Trait",
    y = "Percentage of genomes (%)",
    title = "Distribution of PGP functional traits in SABC genomes"
  ) +
  scale_fill_manual(
    values = c(
      "Phosphorus" = "#2171b5",
      "Indole" = "#31a354",
      "Glutathione" = "#e6550d",
      "Glycine betaine" = "#756bb1",
      "Trehalose" = "#636363",
      "PQQ" = "#bd0026",
      "Siderophore" = "#3182bd",
      "Cytokinin" = "#31a354"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave("result/PGP_function_bacteria_barplot.pdf", width = 8, height = 5)
colnames(df)
library(tidyverse)

# genome_func_export：第一列是ID，后面全部是PGP 0/1
df_long <- genome_func_export %>%
  rename(genomeID = ID) %>%   # 把ID列改名叫genomeID
  pivot_longer(
    cols = -genomeID,
    names_to = "PGP_trait",
    values_to = "presence"
  )
df_A <- df_long %>%
  group_by(genomeID) %>%
  summarise(PGP_richness = sum(presence), .groups = "drop") %>%
  count(PGP_richness, name = "Genome_number")
pA <- ggplot(df_A, aes(x = PGP_richness, y = Genome_number)) +
  geom_col(fill = "#105830", width = 0.7) +
  geom_text(aes(label = Genome_number), vjust = -0.3, size =3.5) +
  labs(x = "PGP richness", y = "Genome number") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank() # 移除全部横竖网格线
  )
pA
ggsave("result/sFig5a.tiff",pA, width = 15, height = 10, dpi = 300)
ggsave("result/sFig5a.jpg",pA, width = 15, height = 10, dpi = 300)
ggsave("result/sFig5a.pdf",pA, width = 15, height = 10, dpi = 300)
library(tidyverse)
library(ComplexUpset)
df_pgp <- read_csv("result/PGP_function.bacteria_01_matrix.csv")   # 8个PGP TRUE/FALSE
df_tax <- read_csv("data/SGB_taxonomy.csv") # User_genome, Family, Phylum等

# 定义颜色向量
phylum_color <- c(
  "Acidobacteriota"="#b6d7a8",
  "Actinomycetota"="#b4a79b",
  "Bacillota"="#f4cccc",
  "Bacteroidota"="#6fa8dc",
  "Chloroflexota"="#e5efd9",
  "Cyanobacteriota"="#c27ba0",
  "Desulfobacterota"="#f8cbad",
  "Gemmatimonadota"="#76a5af",
  "Myxococcota"="#e6b89c",
  "Others"="#dddddd",
  "Planctomycetota"="#ea9999",
  "Pseudomonadota"="#d0c9c0",
  "Verrucomicrobiota"="#996644"
)
colnames(df_tax)
colnames(df_pgp)

# ==========合并PGP矩阵与分类信息，不再生成Nutrient/Growth/Stress==========
df_merge <- df_pgp %>%
  left_join(df_tax, by = "User_genome") %>%
  mutate(
    `Phylum (gtdb)` = str_trim(`Phylum (gtdb)`),
    `Phylum (gtdb)` = replace_na(`Phylum (gtdb)`, "Others"),
    `Phylum (gtdb)` = ifelse(!`Phylum (gtdb)` %in% names(phylum_color), "Others", `Phylum (gtdb)`)
  )

pgp_cols <- c("Phosphorus","Siderophore","Indole","Cytokinin","Glutathione","Trehalose","Glycine_betaine")

pgp_cols <- c(
  "Phosphorus",
  "Siderophore",
  "Indole",
  "Cytokinin",
  "Glutathione",
  "Trehalose",
  "Glycine_betaine"
)

p_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Phosphorus = "Phosphorus",
      Siderophore = "Siderophore",
      Indole = "IAA",
      Cytokinin = "CK",
      Glutathione = "Glutathione",
      Trehalose = "Trehalose",
      Glycine_betaine = "Glycine betaine"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

p_upset
ggsave("result/sFig5b.tiff",p_upset, width = 15, height = 10, dpi = 300)
ggsave("result/sFig5b.jpg",p_upset, width = 15, height = 10, dpi = 300)
ggsave("result/sFig5b.pdf",p_upset, width = 15, height = 10, dpi = 300)
colnames(df_merge)

pgp_cols <- c("Phosphorus","Siderophore")

pgp_cols <- c(
  "Phosphorus",
  "Siderophore"
)

pC_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Phosphorus = "Phosphorus",
      Siderophore = "Siderophore"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

pC_upset
ggsave("result/sFig5c.tiff",pC_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig5c.jpg",pC_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig5c.pdf",pC_upset, width = 10, height = 10, dpi = 300)

pgp_cols <- c("Indole","Cytokinin")

pgp_cols <- c(
  "Indole",
  "Cytokinin"
)

pD_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Indole = "IAA",
      Cytokinin = "CK"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

pD_upset
ggsave("result/sFig5d.tiff",pD_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig5d.jpg",pD_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig5d.pdf",pD_upset, width = 10, height = 10, dpi = 300)

pgp_cols <- c("Glutathione","Trehalose","Glycine_betaine")

pgp_cols <- c(
  "Glutathione",
  "Trehalose",
  "Glycine_betaine"
)

pE_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Glutathione = "Glutathione",
      Trehalose = "Trehalose",
      Glycine_betaine = "Glycine betaine"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

pE_upset
ggsave("result/sFig5e.tiff",pE_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig5e.jpg",pE_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig5e.pdf",pE_upset, width = 10, height = 10, dpi = 300)

# 第一行：pA | p_upset；第二行：pC_upset | pD_upset | pE_upset
sFig5 <- (pA | p_upset) / (pC_upset | pD_upset | pE_upset) +
  plot_annotation(tag_levels = "a") +
  plot_layout(
    guides = "collect",
    heights = c(1, 1),    # 上下两行高度比例
    widths = c(1,1,1,1,1) # 5个图各自宽度权重
  )
ggsave("result/sFig5.tiff",sFig5, width = 20, height = 15, dpi = 300)
ggsave("result/sFig5.jpg",sFig5, width = 20, height = 15, dpi = 300)
ggsave("result/sFig5.pdf",sFig5, width = 20, height = 15, dpi = 300)

library(tidyverse)

# ====================== 1. 读取数据+基础配置 ======================
anno_df <- read.csv("data/Archaea.Eggnog.csv", stringsAsFactors = F)
meta <- read.csv("data/SGB_taxonomy.csv", stringsAsFactors = F)
# ID列就是基因组/bin名称！
genome_col <- "ID"  

# 8个目标功能定义
assign_function <- function(desc){
  if(is.na(desc)) return(NA)
  desc_low <- tolower(desc)
  if(str_detect(desc_low, "phosphorus|phosphate|phosphonate")) return("Phosphorus")
  if(str_detect(desc_low, "siderophore")) return("Siderophore")
  if(str_detect(desc_low, "indole")) return("Indole")
  if(str_detect(desc_low, "cytokinin")) return("Cytokinin")
  if(str_detect(desc_low, "glutathione")) return("Glutathione")
  if(str_detect(desc_low, "trehalose")) return("Trehalose")
  if(str_detect(desc_low, "glycine betaine|glycine-betaine|betaine")) return("Glycine betaine")
  if(str_detect(desc_low, "pqq")) return("PQQ")
  return(NA)
}

# ====================== 2. 清洗+整理功能-KO映射表 ======================
anno_clean <- anno_df %>%
  mutate(Function = map_chr(Description, assign_function)) %>%
  filter(!is.na(Function), KEGG_ko != "-", !is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_remove(KEGG_ko, "^ko:")) %>%
  distinct(Function, KEGG_ko, .data[[genome_col]])
func_ko_set <- anno_clean %>%
  distinct(Function, KEGG_ko) %>%
  group_by(Function) %>%
  summarise(all_ko = list(KEGG_ko), ko_num = n()) %>%
  ungroup()

func_ko_set
# ====================== 替换这段：逐基因组判定PGP功能存在性 ======================
min_count <- 1  # 至少检出1个该功能KO就算存在

genome_func_matrix <- anno_clean %>%
  group_by(ID) %>%
  summarise(genome_ko = list(unique(KEGG_ko))) %>%
  rowwise() %>%
  mutate(
    Phosphorus = {
      ko_target <- func_ko_set %>% filter(Function == "Phosphorus") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Siderophore = {
      ko_target <- func_ko_set %>% filter(Function == "Siderophore") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Indole = {
      ko_target <- func_ko_set %>% filter(Function == "Indole") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Cytokinin = {
      ko_target <- func_ko_set %>% filter(Function == "Cytokinin") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Glutathione = {
      ko_target <- func_ko_set %>% filter(Function == "Glutathione") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Trehalose = {
      ko_target <- func_ko_set %>% filter(Function == "Trehalose") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    Glycine_betaine = {
      ko_target <- func_ko_set %>% filter(Function == "Glycine betaine") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    },
    PQQ = {
      ko_target <- func_ko_set %>% filter(Function == "PQQ") %>% pull(all_ko) %>% unlist()
      sum(ko_target %in% genome_ko) >= min_count
    }
  ) %>%
  ungroup()

# 保存0/1矩阵（去掉list列genome_ko！！）
genome_func_export <- genome_func_matrix %>% select(-genome_ko)

write.csv(genome_func_export, "result/PGP_function.archaea_01_matrix.csv", row.names = F)

# ====================== 统计检出比例+绘图（直接跑） ======================
total_genome <- nrow(genome_func_export)

func_summary <- genome_func_export %>%
  summarise(across(!ID, sum)) %>%  # !ID 表示除了ID以外所有列
  pivot_longer(everything(), names_to = "Function", values_to = "Present_genome") %>%
  mutate(
    Ratio = Present_genome / total_genome * 100,
    Function = str_replace(Function, "_", " ")
  ) %>%
  arrange(desc(Ratio))

print(func_summary)
write.csv(func_summary, "result/PGP_function_archaea_ratio.csv", row.names = F)

# 绘图
ggplot(func_summary, aes(x = reorder(Function, Ratio), y = Ratio, fill = Function)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(Ratio,1), "%")), hjust = -0.2, size = 4) +
  coord_flip() +
  labs(
    x = "PGP Functional Trait",
    y = "Percentage of genomes (%)",
    title = "Distribution of PGP functional traits in SABC genomes"
  ) +
  scale_fill_manual(
    values = c(
      "Phosphorus" = "#2171b5",
      "Indole" = "#31a354",
      "Glutathione" = "#e6550d",
      "Glycine betaine" = "#756bb1",
      "Trehalose" = "#636363",
      "PQQ" = "#bd0026",
      "Siderophore" = "#3182bd",
      "Cytokinin" = "#31a354"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave("result/PGP_function_archaea_barplot.pdf", width = 8, height = 5)
colnames(df)
library(tidyverse)

# genome_func_export：第一列是ID，后面全部是PGP 0/1
df_long <- genome_func_export %>%
  rename(genomeID = ID) %>%   # 把ID列改名叫genomeID
  pivot_longer(
    cols = -genomeID,
    names_to = "PGP_trait",
    values_to = "presence"
  )
df_A <- df_long %>%
  group_by(genomeID) %>%
  summarise(PGP_richness = sum(presence), .groups = "drop") %>%
  count(PGP_richness, name = "Genome_number")
pA <- ggplot(df_A, aes(x = PGP_richness, y = Genome_number)) +
  geom_col(fill = "#aaabc2", width = 0.7) +
  geom_text(aes(label = Genome_number), vjust = -0.3, size =3.5) +
  labs(x = "PGP richness", y = "Genome number") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank() # 移除全部横竖网格线
  )
pA
ggsave("result/sFig6a.tiff",pA, width = 15, height = 10, dpi = 300)
ggsave("result/sFig6a.jpg",pA, width = 15, height = 10, dpi = 300)
ggsave("result/sFig6a.pdf",pA, width = 15, height = 10, dpi = 300)
library(tidyverse)
library(ComplexUpset)
df_pgp <- read_csv("result/PGP_function.archaea_01_matrix.csv")   # 8个PGP TRUE/FALSE
df_tax <- read_csv("data/SGB_taxonomy.csv") # User_genome, Family, Phylum等

phylum_color <- c(
  "Halobacteriota" = "#4a5f76",
  "Thermoplasmatota" = "#4c6b57",
  "Thermoproteota" = "#d9b38c",
  "Nanoarchaeota"="#b0cfc7",
  "Nanohaloarchaeota"="#aaabc2",
  "SpSt-1190"="#cae7df"
)
colnames(df_tax)
colnames(df_pgp)

# ==========合并PGP矩阵与分类信息，不再生成Nutrient/Growth/Stress==========
df_merge <- df_pgp %>%
  left_join(df_tax, by = "User_genome") %>%
  mutate(
    `Phylum (gtdb)` = str_trim(`Phylum (gtdb)`),
    `Phylum (gtdb)` = replace_na(`Phylum (gtdb)`, "Others"),
    `Phylum (gtdb)` = ifelse(!`Phylum (gtdb)` %in% names(phylum_color), "Others", `Phylum (gtdb)`)
  )
# 统计每个门，包括NA
df_merge %>% count(`Phylum (gtdb)`)
df_merge %>%
  filter(`Phylum (gtdb)` == "Others") %>%
  select(User_genome, `Phylum (gtdb)`, `Class (gtdb)`, `Order (gtdb)`)
pgp_cols <- c("Phosphorus","Siderophore","Indole","Cytokinin","Glutathione","Trehalose","Glycine_betaine")

pgp_cols <- c(
  "Phosphorus",
  "Siderophore",
  "Indole",
  "Cytokinin",
  "Glutathione",
  "Trehalose",
  "Glycine_betaine"
)

p_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Phosphorus = "Phosphorus",
      Siderophore = "Siderophore",
      Indole = "IAA",
      Cytokinin = "CK",
      Glutathione = "Glutathione",
      Trehalose = "Trehalose",
      Glycine_betaine = "Glycine betaine"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

p_upset
ggsave("result/sFig6b.tiff",p_upset, width = 15, height = 10, dpi = 300)
ggsave("result/sFig6b.jpg",p_upset, width = 15, height = 10, dpi = 300)
ggsave("result/sFig6b.pdf",p_upset, width = 15, height = 10, dpi = 300)
colnames(df_merge)

pgp_cols <- c("Phosphorus","Siderophore")

pgp_cols <- c(
  "Phosphorus",
  "Siderophore"
)

pC_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Phosphorus = "Phosphorus",
      Siderophore = "Siderophore"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

pC_upset
ggsave("result/sFig6c.tiff",pC_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig6c.jpg",pC_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig6c.pdf",pC_upset, width = 10, height = 10, dpi = 300)

pgp_cols <- c("Indole","Cytokinin")

pgp_cols <- c(
  "Indole",
  "Cytokinin"
)

pD_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Indole = "IAA",
      Cytokinin = "CK"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

pD_upset
ggsave("result/sFig6d.tiff",pD_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig6d.jpg",pD_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig6d.pdf",pD_upset, width = 10, height = 10, dpi = 300)

pgp_cols <- c("Glutathione","Trehalose","Glycine_betaine")

pgp_cols <- c(
  "Glutathione",
  "Trehalose",
  "Glycine_betaine"
)

pE_upset <- upset(
  df_merge,
  intersect = pgp_cols,
  min_size = 1,
  width_ratio = 0.25,
  sort_sets = FALSE,
  sort_intersections = FALSE,
  # 修改左侧集合的显示标签，IAA / CK
  set_sizes = upset_set_size() +
    scale_x_discrete(labels = c(
      Glutathione = "Glutathione",
      Trehalose = "Trehalose",
      Glycine_betaine = "Glycine betaine"
    )),
  base_annotations = list(
    'Intersection size' = intersection_size(
      mapping = aes(fill = `Phylum (gtdb)`),
      position = position_stack()
    ) + scale_fill_manual(values = phylum_color, name = "Phylum (gtdb)")
  ),
  themes = upset_modify_themes(list(
    'intersections_matrix' = theme(panel.grid = element_blank()),
    'Intersection size' = theme(panel.grid = element_blank()),
    'overall_sizes' = theme(panel.grid = element_blank())
  ))
) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.7, "line"),
    legend.title = element_text(size = 11, face = "bold")
  )

pE_upset
ggsave("result/sFig6e.tiff",pE_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig6e.jpg",pE_upset, width = 10, height = 10, dpi = 300)
ggsave("result/sFig6e.pdf",pE_upset, width = 10, height = 10, dpi = 300)

# 第一行：pA | p_upset；第二行：pC_upset | pD_upset | pE_upset
sFig6 <- (pA | p_upset) / (pC_upset | pD_upset | pE_upset) +
  plot_annotation(tag_levels = "a") +
  plot_layout(
    guides = "collect",
    heights = c(1, 1),    # 上下两行高度比例
    widths = c(1,1,1,1,1) # 5个图各自宽度权重
  )
sFig6
ggsave("result/sFig6.tiff",sFig6, width = 20, height = 15, dpi = 300)
ggsave("result/sFig6.jpg",sFig6, width = 20, height = 15, dpi = 300)
ggsave("result/sFig6.pdf",sFig6, width = 20, height = 15, dpi = 300)

library(ape)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)

tree <- read.tree("data/gtdbtk_iqtree2_result.treefile") 
metadata <- read.csv("data/SGB.info.csv"
                     )

# ============ 第一步：预处理metadata，自动归类Others（放在最前面！！） ============
# 定义你想要单独上色的门列表
valid_phylum <- c(
  "Acidobacteriota",
  "Actinomycetota",
  "Bacillota",
  "Bacteroidota",
  "Chloroflexota",
  "Cyanobacteriota",
  "Desulfobacterota",
  "Gemmatimonadota",
  "Myxococcota",
  "Planctomycetota",
  "Pseudomonadota",
  "Verrucomicrobiota"
)
colnames(metadata)
metadata <- metadata %>%
  mutate(
    Phylum..gtdb. = case_when(
      Phylum..gtdb. %in% valid_phylum ~ Phylum..gtdb.,
      TRUE ~ "Others" # 剩下全部自动转为Others，包括NA
    ),
    log10_bgc = log10(`BGC_count` + 1)
  )

# ============ 第二步：更新颜色向量，给Others指定颜色（比如#888888深灰，可自行替换） ============
tax_color <- c(
  "Acidobacteriota"="#b6d7a8",
  "Actinomycetota"="#b4a79b",
  "Bacillota"="#f4cccc",
  "Bacteroidota"="#6fa8dc",
  "Chloroflexota"="#e5efd9",
  "Cyanobacteriota"="#c27ba0",
  "Desulfobacterota"="#f8cbad",
  "Gemmatimonadota"="#76a5af",
  "Myxococcota"="#e6b89c",
  "Planctomycetota"="#ea9999",
  "Pseudomonadota"="#d0c9c0",
  "Verrucomicrobiota"="#996644",
  "Others"="#dddddd" # 👈 修改这里换Others颜色，例如#636363深灰，#bdbdbd浅灰
)
colnames(metadata)
# ============ 第三步：原来的树代码不变，直接运行 ============
tree.plot <- ggtree(tree, layout = "fan", right = TRUE, open.angle = 180, size = 0.1) %<+% metadata +
  geom_aline(aes(color = Phylum..gtdb.), linetype = "solid", size = 0.2, show.legend = TRUE) +
  scale_color_manual(values = tax_color, name = "phylum") + # color用tax_color，Others自动匹配
  guides(
    colour = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
    )
  ) +
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.5, "mm")
  )
tree.plot
library(ggtreeExtra)
# ============ 第四步：追加外圈柱子 ============
p_a <- tree.plot +
  geom_fruit(
    geom = geom_col,
    mapping = aes(y = label, x = log10_bgc, fill = Phylum..gtdb.),
    width = 0.65,
    pwidth = 0.22,
    offset = 0.02,
    orientation = "y"
  ) +
  scale_fill_manual(values = tax_color, name = "Phylum", drop = FALSE) +
  guides(
    color = "none",           # 隐藏内侧aline的图例
    fill = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12)
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.2, "mm")
  )
# fill复用同一个tax_color
p_a
ggsave("result/Fig3a.1750.tiff",p_a, width = 10, height = 6, dpi = 300)
ggsave("result/Fig3a.1750.jpg",p_a, width = 10, height = 6, dpi = 300)
ggsave("result/Fig3a.1750.pdf",p_a, width = 10, height = 6, dpi = 300)

library(tidyverse)

df <- read.csv("data/BGC-1750.info.csv") 
keep_phyla <- c(
  "Acidobacteriota",
  "Pseudomonadota",
  "Actinomycetota",
  "Bacteroidota",
  "Myxococcota",
  "Planctomycetota",
  "Gemmatimonadota",
  "Chloroflexota",
  "KSB1",
  "Bacillota"
)

# ⚠️ 全部替换为普通短横杠 "-"，删掉长破折号
bgc_type_col <- c(
  "phosphonate"       = "#ecd8b8",
  "ladderane"         = "#82a7c6",
  "resorcinol"        = "#a7b997",
  "Others"            = "#d2b4b4",
  "NI-siderophore"    = "#b8c8b0",
  "PKS-NRP_Hybrids"   = "#c9b29c",
  "arylpolyene"       = "#708db1",
  "PKSI"              = "#d8c7b0",
  "PKSother"          = "#e7d8c3",
  "RiPPs"             = "#c8c392",
  "Terpene"           = "#9faa92",
  "NRPS"              = "#79a778"
)

genome_info <- df %>%
  distinct(Bacteria_ID, Bacteria_phylum) %>%
  rename(genome_id = Bacteria_ID, phylum = Bacteria_phylum)

bgc_type_long <- df %>%
  count(Bacteria_ID, BiG.SCAPE.class, name = "count") %>%
  rename(genome_id = Bacteria_ID, bgc_type = BiG.SCAPE.class) %>%
  mutate(bgc_type = factor(bgc_type, levels = names(bgc_type_col)))

bgc_phylum <- bgc_type_long %>%
  left_join(genome_info %>% select(genome_id, phylum), by="genome_id") %>%
  mutate(
    phylum = if_else(phylum %in% keep_phyla, phylum, "Others")
  ) %>%
  group_by(phylum, bgc_type) %>%
  summarise(total = sum(count), .groups="drop") %>%
  group_by(phylum) %>%
  mutate(percent = total / sum(total)*100) %>%
  ungroup()

# ⚠️ Y轴修正：rev()反转，BGC数量最大的门出现在Y轴【最上方】
phy_ord <- bgc_phylum %>%
  group_by(phylum) %>%
  summarise(s = sum(total), .groups="drop") %>%
  arrange(desc(s)) %>%
  pull(phylum) %>%
  rev()   # 增加rev()，解决Y轴上下颠倒

bgc_phylum$phylum <- factor(bgc_phylum$phylum, levels = phy_ord)

# Y轴标签：Phylum(总数)
phy_label_df <- bgc_phylum %>%
  distinct(phylum) %>%
  left_join(
    bgc_phylum %>% group_by(phylum) %>% summarise(tot = sum(total)),
    by = "phylum"
  ) %>%
  mutate(phy_label = paste0(phylum,"(",tot,")"))
phy_label_map <- set_names(phy_label_df$phy_label, phy_label_df$phylum)

p_b <- ggplot(bgc_phylum, aes(x = percent, y = phylum)) +
  geom_col(
    aes(fill = bgc_type),
    position = position_stack(),
    width = 0.8, colour = "white", linewidth = 0.1
  ) +
  scale_fill_manual(
    values = bgc_type_col,
    name = "BGC type",
    drop = FALSE,
    limits = names(bgc_type_col)
  ) +
  scale_y_discrete(labels = phy_label_map) +
  labs(
    x = "Percentage of BGCs within Phylum (%)",
    y = "",
    title = "b"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

p_b
ggsave("result/Fig3b.1750.tiff",p_b, width = 9, height = 10, dpi = 300)
ggsave("result/Fig3b.1750.jpg",p_b, width = 9, height = 10, dpi = 300)
ggsave("result/Fig3b.1750.pdf",p_b, width = 9, height = 10, dpi = 300)

library(tidyverse)

keep_phyla <- c(
  "Acidobacteriota",
  "Pseudomonadota",
  "Actinomycetota",
  "Bacteroidota",
  "Chloroflexota",
  "Planctomycetota",
  "Gemmatimonadota",
  "Desulfobacterota",
  "Myxococcota",
  "Bacillota",
  "Others"
)

phylum_color <- c(
  "Acidobacteriota"     = "#b8d8b3",
  "Pseudomonadota"      = "#d0c9c0",
  "Actinomycetota"      = "#b4a79b",
  "Bacteroidota"        = "#7688a5",
  "Chloroflexota"       = "#e2e6c8",
  "Planctomycetota"     = "#e8b3b8",
  "Gemmatimonadota"     = "#8cb4d8",
  "Desulfobacterota"    = "#f7ddbc",
  "Myxococcota"         = "#b08c6e",
  "Bacillota"           = "#76a5af",
  "Others"              = "#dddddd"
)

#==================== 数据处理 ====================
genome_info <- df %>%
  distinct(Bacteria_ID, Bacteria_phylum) %>%
  rename(genome_id = Bacteria_ID, phylum = Bacteria_phylum)

# 统计：每个BGC type × phylum 的count
d_data <- df %>%
  rename(genome_id = Bacteria_ID, bgc_type = BiG.SCAPE.class) %>%
  left_join(genome_info, by = "genome_id") %>%
  mutate(
    phylum = if_else(phylum %in% keep_phyla, phylum, "Others"),
    # =========核心改动：phylum因子顺序 = 图例从上到下顺序，堆叠从左到右固定！=========
    phylum = factor(phylum, levels = names(phylum_color))
  ) %>%
  count(bgc_type, phylum, name = "cnt") %>%
  group_by(bgc_type) %>%
  mutate(
    total_bgc = sum(cnt),
    percent = cnt / sum(cnt) *100
  ) %>%
  ungroup()

# 构造Y轴标签：BGCtype(总数)
y_label_df <- d_data %>%
  distinct(bgc_type, total_bgc) %>%
  mutate(y_label = paste0(bgc_type,"(",total_bgc,")"))
y_map <- set_names(y_label_df$y_label, y_label_df$bgc_type)

# Y轴顺序：总BGC数最大NRPS放图【最上方】，最小phosphonate放最下面，和参考d图对齐
y_level_order <- y_label_df %>%
  arrange(desc(total_bgc)) %>%
  pull(bgc_type) %>%
  rev()

d_data$bgc_type <- factor(d_data$bgc_type, levels = y_level_order)

#==================== 绘图d图 ====================
p_d <- ggplot(d_data, aes(x = percent, y = bgc_type)) +
  # 删掉 group = stack_group，fill直接用phylum因子，全局固定堆叠顺序
  geom_col(
    aes(fill = phylum),
    position = position_stack(),
    width = 0.75, colour = "white", linewidth =0.1
  ) +
  scale_fill_manual(
    values = phylum_color,
    name = "Phylum",
    drop = FALSE,
    limits = names(phylum_color) # 图例顺序锁定，和phylum_color顺序一致
  ) +
  scale_y_discrete(labels = y_map) +
  scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,75,100)) +
  labs(
    x = "Number of Phylum(%)",
    y = "",
    title = "d"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size=8),
    plot.title = element_text(hjust=0, face="bold")
  )

p_d
ggsave("result/Fig3d.1750.tiff",p_d, width = 9, height = 10, dpi = 300)
ggsave("result/Fig3d.1750.jpg",p_d, width = 9, height = 10, dpi = 300)
ggsave("result/Fig3d.1750.pdf",p_d, width = 9, height = 10, dpi = 300)

df <- read.csv("data/Pseudo.BGC.csv")
library(tidyverse)
library(patchwork)

genome_bgc_richness <- df %>%
  group_by(Bacteria_ID, Family) %>%
  summarise(bgc_richness = n(), .groups = "drop")

family_genome_cnt <- genome_bgc_richness %>%
  group_by(Family) %>%
  summarise(genome_n = n_distinct(Bacteria_ID), .groups = "drop")

# --------过滤，保留基因组>=2的Family，避免单点箱线--------
family_genome_cnt_filtered <- family_genome_cnt %>% filter(genome_n >=2)
keep_family <- family_genome_cnt_filtered$Family

genome_bgc_richness_filtered <- genome_bgc_richness %>%
  filter(Family %in% keep_family)

# ✅【关键：只在这里定义一次factor顺序，两张图共用】
family_order <- family_genome_cnt_filtered %>%
  arrange(desc(genome_n)) %>%
  pull(Family)

# 全部数据集统一赋值factor，只定义一次！！
genome_bgc_richness_filtered$Family <- factor(genome_bgc_richness_filtered$Family, levels = family_order)
family_genome_cnt_filtered$Family <- factor(family_genome_cnt_filtered$Family, levels = family_order)

# 左图
p_left <- ggplot(family_genome_cnt_filtered, aes(x = genome_n, y = Family)) +
  geom_bar(stat = "identity", fill = "#309930") +
  scale_x_reverse() +
  labs(x = "Genome number", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p_left
# 右图
p_right <- ggplot(genome_bgc_richness_filtered, aes(x = bgc_richness, y = Family)) +
  geom_boxplot(fill = "#70c870", outlier.shape = NA) +
  labs(x = "BGC richness", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_right
# ✅关键参数：align = "y"，强制Y轴对齐！！
p_c <- p_left + p_right +
  plot_layout(widths = c(1,1.2)) 
p_c
ggsave("result/Fig3c.Pseu.tiff",p_c, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3c.Pseu.jpg",p_c, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3c.Pseu.pdf",p_c, width = 10, height = 10, dpi = 300)

df <- read.csv("data/Acido.BGC.csv")
library(tidyverse)
library(patchwork)

genome_bgc_richness <- df %>%
  group_by(Bacteria_ID, Family) %>%
  summarise(bgc_richness = n(), .groups = "drop")

family_genome_cnt <- genome_bgc_richness %>%
  group_by(Family) %>%
  summarise(genome_n = n_distinct(Bacteria_ID), .groups = "drop")

# --------过滤，保留基因组>=2的Family，避免单点箱线--------
family_genome_cnt_filtered <- family_genome_cnt %>% filter(genome_n >=10)
keep_family <- family_genome_cnt_filtered$Family

genome_bgc_richness_filtered <- genome_bgc_richness %>%
  filter(Family %in% keep_family)

# ✅【关键：只在这里定义一次factor顺序，两张图共用】
family_order <- family_genome_cnt_filtered %>%
  arrange(desc(genome_n)) %>%
  pull(Family)

# 全部数据集统一赋值factor，只定义一次！！
genome_bgc_richness_filtered$Family <- factor(genome_bgc_richness_filtered$Family, levels = family_order)
family_genome_cnt_filtered$Family <- factor(family_genome_cnt_filtered$Family, levels = family_order)

# 左图
p_left <- ggplot(family_genome_cnt_filtered, aes(x = genome_n, y = Family)) +
  geom_bar(stat = "identity", fill = "#309930") +
  scale_x_reverse() +
  labs(x = "Genome number", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p_left
# 右图
p_right <- ggplot(genome_bgc_richness_filtered, aes(x = bgc_richness, y = Family)) +
  geom_boxplot(fill = "#70c870", outlier.shape = NA) +
  labs(x = "BGC richness", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_right
# ✅关键参数：align = "y"，强制Y轴对齐！！
p_out <- p_left + p_right +
  plot_layout(widths = c(1,1.2)) 
p_out
ggsave("result/Fig3c.Acido.tiff",p_out, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3c.Acido.jpg",p_out, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3c.Acido.pdf",p_out, width = 10, height = 10, dpi = 300)

df <- read.csv("data/300.BGC.Acido.csv")
library(tidyverse)
library(patchwork)

genome_bgc_richness <- df %>%
  group_by(Bacteria_ID, Family) %>%
  summarise(bgc_richness = n(), .groups = "drop")

family_genome_cnt <- genome_bgc_richness %>%
  group_by(Family) %>%
  summarise(genome_n = n_distinct(Bacteria_ID), .groups = "drop")

# --------过滤，保留基因组>=2的Family，避免单点箱线--------
family_genome_cnt_filtered <- family_genome_cnt %>% filter(genome_n >=2)
keep_family <- family_genome_cnt_filtered$Family

genome_bgc_richness_filtered <- genome_bgc_richness %>%
  filter(Family %in% keep_family)

# ✅【关键：只在这里定义一次factor顺序，两张图共用】
family_order <- family_genome_cnt_filtered %>%
  arrange(desc(genome_n)) %>%
  pull(Family)

# 全部数据集统一赋值factor，只定义一次！！
genome_bgc_richness_filtered$Family <- factor(genome_bgc_richness_filtered$Family, levels = family_order)
family_genome_cnt_filtered$Family <- factor(family_genome_cnt_filtered$Family, levels = family_order)

# 左图
p_left <- ggplot(family_genome_cnt_filtered, aes(x = genome_n, y = Family)) +
  geom_bar(stat = "identity", fill = "#309930") +
  scale_x_reverse() +
  labs(x = "Genome number", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p_left
# 右图
p_right <- ggplot(genome_bgc_richness_filtered, aes(x = bgc_richness, y = Family)) +
  geom_boxplot(fill = "#70c870", outlier.shape = NA) +
  labs(x = "BGC richness", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_right
# ✅关键参数：align = "y"，强制Y轴对齐！！
p_out <- p_left + p_right +
  plot_layout(widths = c(1,1.2)) 
p_out
ggsave("result/Fig3c.300.Acido.tiff",p_out, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3c.300.Acido.jpg",p_out, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3c.300.Acido.pdf",p_out, width = 10, height = 10, dpi = 300)

#读入数据
df <- read.csv("data/SGB.info.csv")

library(tidyverse)
library(patchwork)

# 门水平配色
tax_color <- c(
  "Acidobacteriota"="#b6d7a8",
  "Actinomycetota"="#b4a79b",
  "Bacillota"="#f4cccc",
  "Bacteroidota"="#6fa8dc",
  "Chloroflexota"="#e5dfff",
  "Cyanobacteriota"="#c27ba0",
  "Desulfobacterota"="#f8cbad",
  "Gemmatimonadota"="#76a5af",
  "Myxococcota"="#e6b89c",
  "Planctomycetota"="#ea9999",
  "Pseudomonadota"="#d0c9c0",
  "Verrucomicrobiota"="#996644",
  "Halobacteriota" = "#4a5f76",
  "Thermoplasmatota" = "#4c6b57",
  "Thermoproteota" = "#d9b38c",
  "Nanoarchaeota"="#b0cfc7",
  "Nanohaloarchaeota"="#aaabc2",
  "Others"="#dddddd"
)

# 预处理：不在配色里的门统一归为Others
df <- df %>%
  mutate(
    genome_size_Mb = size / 1e6,
    bgc_percent = (BGC_count / size) * 1e6 * 100,
    Phylum_clean = ifelse(`Phylum..gtdb.` %in% names(tax_color), `Phylum..gtdb.`, "Others")
  ) %>%
  filter(BGC_count > 0)   # 去掉BGC=0的点，可选，不需要就删掉这行

# 上方密度图
p_density <- ggplot(df, aes(x = bgc_percent)) +
  geom_density(fill = "#d8b898", alpha = 0.7) +
  labs(x = NULL, y = "Density") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p_density
install.packages("ggrepel")
library(ggrepel)
library(tidyverse)
top_label <- df %>% slice_max(bgc_percent, n = 5)

p_scatter <- ggplot(df, aes(x = genome_size_Mb, y = bgc_percent, color = Phylum_clean)) +
  geom_point(size = 1.3, alpha = 0.7) +
  #替换原来geom_label → geom_label_repel自动排斥，避免重叠
  geom_label_repel(
    data = top_label,
    aes(label = paste0("Bacteria: ",User_genome,"\nPhylum: ",Phylum_clean,"\nBGCs: ",BGC_count)),
    size =2.2,
    label.size =0.2,
    label.padding = unit(0.2,"lines"),
    max.overlaps =20,   #允许最多多少标签
    box.padding =0.4,   #标签之间空隙
    point.padding =0.2,
    segment.size=0.2    #指向点的引线粗细
  ) +
  scale_color_manual(values = tax_color, name = "Taxonomy") +
  labs(
    x = "Genome size (Mb)",
    y = "BGC encoding percentage (%)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face="bold"),
    legend.key.spacing.y = unit(1.2,"mm")
  )
p_scatter
p_e <- p_density / p_scatter + plot_layout(heights = c(1, 4))
p_e
ggsave("result/Fig3e.tiff",p_e, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3e.jpg",p_e, width = 10, height = 10, dpi = 300)
ggsave("result/Fig3e.pdf",p_e, width = 10, height = 10, dpi = 300)

library(ape)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)

tree <- read.tree("data/gtdbtk.unrooted2.tree") 
metadata <- read.csv("data/SGB.info.csv", 
)

# ============ 第一步：预处理metadata，自动归类Others（放在最前面！！） ============
# 定义你想要单独上色的门列表
valid_phylum <- c(
  "Halobacteriota",
  "Thermoplasmatota",
  "Thermoproteota",
  "Nanoarchaeota",
  "Nanohaloarchaeota",
  "SpSt-1190"
)
colnames(metadata)
metadata <- metadata %>%
  mutate(
    Phylum..gtdb. = case_when(
      Phylum..gtdb. %in% valid_phylum ~ Phylum..gtdb.,
      TRUE ~ "Others" # 剩下全部自动转为Others，包括NA
    ),
    log10_bgc = log10(`BGC_count` + 1)
  )

# ============ 第二步：更新颜色向量，给Others指定颜色（比如#888888深灰，可自行替换） ============
tax_color <- c(
  "Halobacteriota" = "#4a5f76",
  "Thermoplasmatota" = "#4c6b57",
  "Thermoproteota" = "#d9b38c",
  "Nanoarchaeota"="#b0cfc7",
  "Nanohaloarchaeota"="#aaabc2",
  "SpSt-1190"="#cae7df"
)
colnames(metadata)
# ============ 第三步：原来的树代码不变，直接运行 ============
tree.plot <- ggtree(tree, layout = "fan", right = TRUE, open.angle = 180, size = 0.1) %<+% metadata +
  geom_aline(aes(color = Phylum..gtdb.), linetype = "solid", size = 0.2, show.legend = TRUE) +
  scale_color_manual(values = tax_color, name = "phylum") + # color用tax_color，Others自动匹配
  guides(
    colour = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
    )
  ) +
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.5, "mm")
  )
tree.plot
library(ggtreeExtra)
# ============ 第四步：追加外圈柱子 ============
p_a <- tree.plot +
  geom_fruit(
    geom = geom_col,
    mapping = aes(y = label, x = log10_bgc, fill = Phylum..gtdb.),
    width = 0.65,
    pwidth = 0.22,
    offset = 0.02,
    orientation = "y"
  ) +
  scale_fill_manual(values = tax_color, name = "Phylum", drop = FALSE) +
  guides(
    color = "none",           # 隐藏内侧aline的图例
    fill = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12)
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.2, "mm")
  )
# fill复用同一个tax_color
p_a
ggsave("result/sFig9a.1750.tiff",p_a, width = 10, height = 6, dpi = 300)
ggsave("result/sFig9a.1750.jpg",p_a, width = 10, height = 6, dpi = 300)
ggsave("result/sFig9a.1750.pdf",p_a, width = 10, height = 6, dpi = 300)

library(tidyverse)

df <- read.csv("data/BGC-178.info.csv") 
keep_phyla <- c(
  "Halobacteriota",
  "Thermoplasmatota",
  "Thermoproteota",
  "Nanoarchaeota",
  "Nanohaloarchaeota",
  "SpSt-1190")

# ⚠️ 全部替换为普通短横杠 "-"，删掉长破折号
bgc_type_col <- c(
  "phosphonate"       = "#ecd8b8",
  "NI-siderophore"    = "#b8c8b0",
  "PKSI"              = "#d8c7b0",
  "PKSother"          = "#e7d8c3",
  "RiPPs"             = "#c8c392",
  "Terpene"           = "#9faa92",
  "NRPS"              = "#79a778"
)

genome_info <- df %>%
  distinct(Archaea_ID, Archaea_phylum) %>%
  rename(genome_id = Archaea_ID, phylum = Archaea_phylum)

bgc_type_long <- df %>%
  count(Archaea_ID, BiG.SCAPE.class, name = "count") %>%
  rename(genome_id = Archaea_ID, bgc_type = BiG.SCAPE.class) %>%
  mutate(bgc_type = factor(bgc_type, levels = names(bgc_type_col)))

bgc_phylum <- bgc_type_long %>%
  left_join(genome_info %>% select(genome_id, phylum), by="genome_id") %>%
  mutate(
    phylum = if_else(phylum %in% keep_phyla, phylum, "Others")
  ) %>%
  group_by(phylum, bgc_type) %>%
  summarise(total = sum(count), .groups="drop") %>%
  group_by(phylum) %>%
  mutate(percent = total / sum(total)*100) %>%
  ungroup()

# ⚠️ Y轴修正：rev()反转，BGC数量最大的门出现在Y轴【最上方】
phy_ord <- bgc_phylum %>%
  group_by(phylum) %>%
  summarise(s = sum(total), .groups="drop") %>%
  arrange(desc(s)) %>%
  pull(phylum) %>%
  rev()   # 增加rev()，解决Y轴上下颠倒

bgc_phylum$phylum <- factor(bgc_phylum$phylum, levels = phy_ord)

# Y轴标签：Phylum(总数)
phy_label_df <- bgc_phylum %>%
  distinct(phylum) %>%
  left_join(
    bgc_phylum %>% group_by(phylum) %>% summarise(tot = sum(total)),
    by = "phylum"
  ) %>%
  mutate(phy_label = paste0(phylum,"(",tot,")"))
phy_label_map <- set_names(phy_label_df$phy_label, phy_label_df$phylum)

p_b <- ggplot(bgc_phylum, aes(x = percent, y = phylum)) +
  geom_col(
    aes(fill = bgc_type),
    position = position_stack(),
    width = 0.8, colour = "white", linewidth = 0.1
  ) +
  scale_fill_manual(
    values = bgc_type_col,
    name = "BGC type",
    drop = FALSE,
    limits = names(bgc_type_col)
  ) +
  scale_y_discrete(labels = phy_label_map) +
  labs(
    x = "Percentage of BGCs within Phylum (%)",
    y = "",
    title = "b"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8)
  )

p_b
ggsave("result/Fig3b.178.tiff",p_b, width = 9, height = 10, dpi = 300)
ggsave("result/Fig3b.178.jpg",p_b, width = 9, height = 10, dpi = 300)
ggsave("result/Fig3b.178.pdf",p_b, width = 9, height = 10, dpi = 300)

library(tidyverse)

keep_phyla <- c(
  "Halobacteriota",
  "Thermoplasmatota",
  "Thermoproteota",
  "Nanoarchaeota",
  "Nanohaloarchaeota")

phylum_color <- c(
  "Halobacteriota" = "#4a5f76",
  "Thermoplasmatota" = "#4c6b57",
  "Thermoproteota" = "#d9b38c",
  "Nanoarchaeota"="#b0cfc7",
  "Nanohaloarchaeota"="#aaabc2"
)

#==================== 数据处理 ====================
genome_info <- df %>%
  distinct(Archaea_ID, Archaea_phylum) %>%
  rename(genome_id = Archaea_ID, phylum = Archaea_phylum)

# 统计：每个BGC type × phylum 的count
d_data <- df %>%
  rename(genome_id = Archaea_ID, bgc_type = BiG.SCAPE.class) %>%
  left_join(genome_info, by = "genome_id") %>%
  mutate(
    phylum = if_else(phylum %in% keep_phyla, phylum, "Others"),
    # =========核心改动：phylum因子顺序 = 图例从上到下顺序，堆叠从左到右固定！=========
    phylum = factor(phylum, levels = names(phylum_color))
  ) %>%
  count(bgc_type, phylum, name = "cnt") %>%
  group_by(bgc_type) %>%
  mutate(
    total_bgc = sum(cnt),
    percent = cnt / sum(cnt) *100
  ) %>%
  ungroup()

# 构造Y轴标签：BGCtype(总数)
y_label_df <- d_data %>%
  distinct(bgc_type, total_bgc) %>%
  mutate(y_label = paste0(bgc_type,"(",total_bgc,")"))
y_map <- set_names(y_label_df$y_label, y_label_df$bgc_type)

# Y轴顺序：总BGC数最大NRPS放图【最上方】，最小phosphonate放最下面，和参考d图对齐
y_level_order <- y_label_df %>%
  arrange(desc(total_bgc)) %>%
  pull(bgc_type) %>%
  rev()

d_data$bgc_type <- factor(d_data$bgc_type, levels = y_level_order)
anyNA(d_data)
which(apply(d_data,1,anyNA))
#==================== 绘图d图 ====================
p_d <- ggplot(d_data, aes(x = percent, y = bgc_type)) +
  geom_col(
    aes(fill = phylum),
    position = position_stack(),
    width = 0.75, colour = "white", linewidth =0.1
  ) +
  scale_fill_manual(
    values = phylum_color,
    drop = FALSE
  ) +
  scale_y_discrete(labels = y_map) +
  scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,75,100)) +
  labs(
    x = "Number of Phylum(%)",
    y = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size=8),
    plot.title = element_text(hjust=0, face="bold"),
    panel.grid = element_blank()
  )

p_d
ggsave("result/sFig9d.178.tiff",p_d, width = 5, height = 5, dpi = 300)
ggsave("result/sFig9d.178.jpg",p_d, width = 5, height = 5, dpi = 300)
ggsave("result/sFig9d.178.pdf",p_d, width = 5, height = 5, dpi = 300)

df <- read.csv("data/Halo.BGC.csv")
library(tidyverse)
library(patchwork)

genome_bgc_richness <- df %>%
  group_by(Archaea_ID, Family) %>%
  summarise(bgc_richness = n(), .groups = "drop")

family_genome_cnt <- genome_bgc_richness %>%
  group_by(Family) %>%
  summarise(genome_n = n_distinct(Archaea_ID), .groups = "drop")

# --------过滤，保留基因组>=2的Family，避免单点箱线--------
family_genome_cnt_filtered <- family_genome_cnt %>% filter(genome_n >=2)
keep_family <- family_genome_cnt_filtered$Family

genome_bgc_richness_filtered <- genome_bgc_richness %>%
  filter(Family %in% keep_family)

# ✅【关键：只在这里定义一次factor顺序，两张图共用】
family_order <- family_genome_cnt_filtered %>%
  arrange(desc(genome_n)) %>%
  pull(Family)

# 全部数据集统一赋值factor，只定义一次！！
genome_bgc_richness_filtered$Family <- factor(genome_bgc_richness_filtered$Family, levels = family_order)
family_genome_cnt_filtered$Family <- factor(family_genome_cnt_filtered$Family, levels = family_order)

# 左图
p_left <- ggplot(family_genome_cnt_filtered, aes(x = genome_n, y = Family)) +
  geom_bar(stat = "identity", fill = "#309930") +
  scale_x_reverse() +
  labs(x = "Genome number", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p_left
# 右图
p_right <- ggplot(genome_bgc_richness_filtered, aes(x = bgc_richness, y = Family)) +
  geom_boxplot(fill = "#70c870", outlier.shape = NA) +
  labs(x = "BGC richness", y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_right
# ✅关键参数：align = "y"，强制Y轴对齐！！
p_c <- p_left + p_right +
  plot_layout(widths = c(1,1.2)) 
p_c
ggsave("result/sFig9c.Halo.tiff",p_c, width = 6, height = 6, dpi = 300)
ggsave("result/sFig9c.Halo.jpg",p_c, width = 6, height = 6, dpi = 300)
ggsave("result/sFig9c.Halo.pdf",p_c, width = 6, height = 6, dpi = 300)
library(patchwork)

final_fig <- (p_d + p_left + p_right) +
  plot_layout(
    heights = c(1),
    widths = c(0.5,0.2,0.2),
    guides = "collect"  # 收集所有图例，放到整体右边
  )

final_fig
ggsave("result/sFig9.tiff",final_fig, width = 12, height = 12, dpi = 300)
ggsave("result/sFig9.jpg",final_fig, width = 8, height = 12, dpi = 300)
ggsave("result/sFig9.pdf",final_fig, width = 8, height = 12, dpi = 300)

library(tidyverse)
library(plotrix)

#==================== 1.花瓣绘图函数（带极值压缩，消除超长花瓣） ====================
flower_plot <- function(sample_names,
                        accessory_counts,
                        display_labels,
                        core_otu,
                        start = 90,
                        inner_radius = 110,
                        extend_max = 130,
                        lab.cex = 0.52){
  
  n <- length(sample_names)
  deg_step <- 360 / n
  
  old_par <- par(pty = "s", mar = c(1,1,2,1))
  on.exit(par(old_par))
  
  acc <- as.numeric(accessory_counts)
  q90 <- quantile(acc, 0.90, na.rm = TRUE)
  acc_clip <- pmin(acc, q90)
  acc_clip <- pmax(acc_clip, 0)
  max_clip <- max(acc_clip)
  extend_r <- (sqrt(acc_clip) / sqrt(max_clip)) * extend_max
  
  total_max_r <- inner_radius + extend_max
  plot(c(-total_max_r, total_max_r),
       c(-total_max_r, total_max_r),
       type = "n", axes = FALSE, xlab = "", ylab = "")
  
  draw.circle(0,0,radius=inner_radius,col="white",border="black",lwd=1.2)
  text(0,0,labels=paste0("Core: ",core_otu),cex=1.1)
  
  pal <- c("#b6d7a8","#b4a79b","#f4cccc","#6fa8dc","#e5efd9","#c27ba0","#f8cbad","#76a5af",
           "#e6b89c","#ea9999","#d0c9c0","#996644","#8dd3c7","#ffffb3","#bebada","#fb8072")
  pal <- c(
    "#c79fbf", "#d7e2c4", "#a5c2bf", "#e6d2ad",
    "#d5b19f", "#b7c3a4", "#c7b1d2", "#d5d5a4",
    "#b3bea4", "#e0c597", "#c5a5af", "#a3c2a4",
    "#d2b3c2", "#c7d6a2", "#b1a3c2", "#dfd2a3"
  )
  
  for(t in seq_along(sample_names)){
    petal_outer_r <- inner_radius + extend_r[t]
    theta1 <- (start + deg_step*(t-1))*pi/180
    theta2 <- (start + deg_step*t)*pi/180
    theta_mid <- (theta1+theta2)/2
    
    x <- c(inner_radius*cos(theta1), petal_outer_r*cos(theta_mid), inner_radius*cos(theta2), inner_radius*cos(theta1))
    y <- c(inner_radius*sin(theta1), petal_outer_r*sin(theta_mid), inner_radius*sin(theta2), inner_radius*sin(theta1))
    polygon(x,y,col=pal[((t-1)%%length(pal))+1],border="gray70",lwd=1)
    
    txt_r <- petal_outer_r *1.06
    text(txt_r*cos(theta_mid), txt_r*sin(theta_mid), labels=display_labels[t], cex=lab.cex, srt = deg_step*(t-1)+start)
  }
}

#==================== 2.批量循环处理8个门 ====================
# 修改这里为你的8个文件名称
#==================== 2.批量循环处理8个门 【修复版】 ====================
file_list <- c(
  "Orthogroups.GeneCount.Pseudomonadota.tsv",
  "Orthogroups.GeneCount.Actinomycetota.tsv",
  "Orthogroups.GeneCount.Bacteroidota.tsv",
  "Orthogroups.GeneCount.Chloroflexota.tsv",
  "Orthogroups.GeneCount.Acidobacteriota.tsv",
  "Orthogroups.GeneCount-Thermoplasmatota.tsv",
  "Orthogroups.GeneCount.Planctomycetota.tsv",
  "Orthogroups.GeneCount.Gemmatimonadota.tsv"
)
#不要用 add_row，初始化空list存储每一行结果
summary_list <- list()

for(f in file_list){
  tax_name <- str_remove(f,"Orthogroups.GeneCount.") %>% str_remove(".tsv")
  message("\n==== Processing: ",tax_name," ====")
  
  df <- read_tsv(f, show_col_types = F)
  
  og_col <- colnames(df)[1]
  genome_cols <- colnames(df)[-1]
  n_genomes <- length(genome_cols)
  
  df_stat <- df %>%
    rowwise() %>%
    mutate(present = sum(c_across(all_of(genome_cols))>0),
           is_core = present == n_genomes) %>%
    ungroup()
  
  core_num <- sum(df_stat$is_core)
  core_og_df <- df_stat %>% filter(is_core)
  accessory_og_df <- df_stat %>% filter(!is_core)
  
  write_csv(core_og_df,paste0(tax_name,"_coreOG.csv"))
  
  count_matrix <- df %>% select(-all_of(og_col)) %>% as.matrix()
  total_per_genome <- colSums(count_matrix>0)
  accessory_per_genome <- pmax(total_per_genome - core_num,0)
  
  #把当前门类结果存成单行tibble，放入list
  current_row <- tibble(
    taxon = tax_name,
    n_genomes = n_genomes,
    core_OG = core_num,
    accessory_OG = nrow(accessory_og_df),
    total_OG = nrow(df_stat)
  )
  summary_list[[tax_name]] <- current_row
  
  pdf(paste0(tax_name,"_pan_flower.pdf"),width=7,height=7)
  flower_plot(
    sample_names = names(accessory_per_genome),
    accessory_counts = accessory_per_genome,
    display_labels = as.integer(total_per_genome),
    core_otu = core_num,
    inner_radius = 110,
    extend_max =130,
    lab.cex =0.2
  )
  title(main=tax_name,cex.main=0.7)
  dev.off()
  
  message(sprintf("%s | Core=%d | Accessory range: %d ~ %d",
                  tax_name,core_num,min(accessory_per_genome),max(accessory_per_genome)))
}

#循环结束，合并全部结果
summary_tb <- bind_rows(summary_list)
write_csv(summary_tb,"pan_genome_summary_all_phyla.csv")
print(summary_tb)

library(tidyverse)
library(patchwork)
library(gridGraphics)

file_list <- c(
  "Orthogroups.GeneCount.Pseudomonadota.tsv",
  "Orthogroups.GeneCount.Actinomycetota.tsv",
  "Orthogroups.GeneCount.Bacteroidota.tsv",
  "Orthogroups.GeneCount.Chloroflexota.tsv",
  "Orthogroups.GeneCount.Acidobacteriota.tsv",
  "Orthogroups.GeneCount-Thermoplasmatota.tsv",
  "Orthogroups.GeneCount.Planctomycetota.tsv",
  "Orthogroups.GeneCount.Gemmatimonadota.tsv"
)

plot_list <- list()

for(f in file_list){
  tax_name <- str_remove(f,"Orthogroups.GeneCount.") %>% str_remove(".tsv")
  message("\n==== Processing: ",tax_name," ====")
  
  df <- read_tsv(f, show_col_types = F)
  og_col <- colnames(df)[1]
  genome_cols <- colnames(df)[-1]
  n_genomes <- length(genome_cols)
  
  df_stat <- df %>%
    rowwise() %>%
    mutate(present = sum(c_across(all_of(genome_cols))>0),
           is_core = present == n_genomes) %>%
    ungroup()
  
  core_num <- sum(df_stat$is_core)
  core_og_df <- df_stat %>% filter(is_core)
  write_csv(core_og_df,paste0(tax_name,"_coreOG.csv"))
  
  count_matrix <- df %>% select(-all_of(og_col)) %>% as.matrix()
  total_per_genome <- colSums(count_matrix>0)
  accessory_per_genome <- pmax(total_per_genome - core_num,0)
  
  # 直接绘图！！不使用recordPlot
  flower_plot(
    sample_names = names(accessory_per_genome),
    accessory_counts = accessory_per_genome,
    display_labels = as.integer(total_per_genome),
    core_otu = core_num,
    inner_radius = 110,
    extend_max = 110,
    lab.cex = 0.2
  )
  title(main=tax_name,cex.main=0.7)
}
dev.off() # 关闭pdf设备，文件写入完成 
library(patchwork)
orthogroups <- read.delim("data/Orthogroups.Chloroflexota.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#e5dfff"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p1 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#e5dfff", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#e5dfff") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Chloroflexota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p1

orthogroups <- read.delim("data/Orthogroups.Acidobacteriota.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#b6d7a8"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p2 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#b6d7a8", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#b6d7a8") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Acidobacteriota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p2

orthogroups <- read.delim("data/Orthogroups.Actinomycetota.tsv", header = TRUE, row.names = 1)
# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#b4a79b"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p3 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#b4a79b", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#b4a79b") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Actinomycetota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p3
orthogroups <- read.delim("data/Orthogroups.Bacteroidota.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#6fa8dc"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p4 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#6fa8dc", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#6fa8dc") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Bacteroidota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p4
orthogroups <- read.delim("data/Orthogroups.Gemmatimonadota.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#76a5af"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p5 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#76a5af", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#76a5af") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Gemmatimonadota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p5

orthogroups <- read.delim("data/Orthogroups.Planctomycetota.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#ea9999"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p6 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#ea9999", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#ea9999") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Planctomycetota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p6

orthogroups <- read.delim("data/Orthogroups.Pseudomonas.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#d0c9c0"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p7 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#d0c9c0", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#d0c9c0") +
  geom_line(aes(y = fit_y), color = "gray", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Pseudomonadota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p7

orthogroups <- read.delim("data/Orthogroups.thermoplas.tsv", header = TRUE, row.names = 1)

# 计算每个家族出现的基因组数量
genome_count <- apply(orthogroups, 1, function(x) sum(x != ""))
n_genomes <- ncol(orthogroups)

# 模拟添加基因组的顺序（随机多次抽样）
n_replicates <- 100
sample_sizes <- 1:n_genomes
pangenome_sizes <- matrix(0, nrow = n_replicates, ncol = n_genomes)

for (rep in 1:n_replicates) {
  order <- sample(n_genomes)
  cumulative <- c()
  current_orthogroups <- c()
  
  for (i in 1:n_genomes) {
    genomes_in <- order[1:i]
    og_in <- apply(orthogroups[, genomes_in, drop = FALSE], 1, function(x) any(x != ""))
    current_orthogroups <- union(current_orthogroups, names(og_in[og_in]))
    cumulative[i] <- length(current_orthogroups)
  }
  
  pangenome_sizes[rep, ] <- cumulative
}

# 计算均值和标准差
mean_pangenome <- apply(pangenome_sizes, 2, mean)
sd_pangenome <- apply(pangenome_sizes, 2, sd)

# 构建数据框
plot_df <- data.frame(
  Genomes = sample_sizes,
  Mean = mean_pangenome,
  SD = sd_pangenome
)

# 拟合 Heaps' law: y = k * x^gamma
heaps_fit <- nls(Mean ~ k * Genomes^gamma, data = plot_df, start = list(k = 100, gamma = 0.5))
k_est <- coef(heaps_fit)["k"]
gamma_est <- coef(heaps_fit)["gamma"]
shadow_color <- "#456b57"
plot_df$fit_y <- predict(heaps_fit)
# 绘图
p8 <- ggplot(plot_df, aes(x = Genomes, y = Mean)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), fill = "#456b57", alpha = 0.4) +
  geom_line(linewidth  = 1.2, color = "#456b57") +
  geom_line(aes(y = fit_y), color = "#dddddd", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.7 * n_genomes, y = max(mean_pangenome), 
           label = paste0("Heaps' law: y = ", round(k_est, 2), " * x^", round(gamma_est, 3)),
           hjust = 0, size = 3) +
  labs(x = "Number of Thermoplasmatota genomes", y = "Pangenome size") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    panel.grid = element_blank()
  )
p8
combined_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, 
                            nrow = 4, ncol = 2, byrow = TRUE)
combined_plot
ggsave("result/sFig10.tiff",combined_plot, width = 12, height = 15, dpi = 300)
ggsave("result/sFig10.jpg",combined_plot, width = 12, height = 15, dpi = 300)
ggsave("result/sFig10.pdf",combined_plot, width = 12, height = 15, dpi = 300)
library(tidyverse)

df <- read_csv("data/8.phylum.cog.csv")
colnames(df)
colnames(df) <- c("Sample", "cog_code", "raw_count", "cog_name", "percent")
colnames(df)
cog_order <- c(
  "Others",
  "Inorganic ion transport and metabolism",
  "Defense mechanisms",
  "Secondary metabolites biosynthesis, transport and catabolism",
  "Signal transduction mechanisms",
  "Cell motility",
  "Cell cycle control, cell division, chromosome partitioning",
  "Cell wall/membrane/envelope biogenesis",
  "Intracellular trafficking, secretion, and vesicular transport",
  "Posttranslational modification, protein turnover, chaperones",
  "Carbohydrate transport and metabolism",
  "Lipid transport and metabolism",
  "Amino acid transport and metabolism",
  "Replication, recombination and repair",
  "Transcription",
  "Coenzyme transport and metabolism",
  "Nucleotide transport and metabolism",
  "Energy production and conversion",
  "Function unknown",
  "Translation, ribosomal structure and biogenesis"
)

df <- df %>%
  mutate(
    Sample = factor(Sample, levels = unique(Sample)),
    COG_full = factor(cog_name, levels = cog_order)
  )
library(ggplot2)

p_stack <- ggplot(df, aes(x = Sample, y = percent, fill = COG_full)) +
  geom_col(position = position_stack()) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Phylum", y = "Proportion (%)", fill = "COG Functional Category") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
p_stack
library(tidyverse)

# 先看全部列名，确认
colnames(df)

# 删除第6列空列，只保留前5列
df <- df[,1:5]

# 现在再重命名5列
colnames(df) <- c("Sample", "cog_code", "raw_count", "cog_name", "percent")

# 再做归一化
df_norm <- df %>%
  group_by(Sample) %>%
  mutate(total = sum(raw_count),
         frac = raw_count / total) %>%
  ungroup()

# 固定COG配色，和你目标参考图匹配
cog_colors <- c(
  "Others" = "#D3D3D3",
  "Translation, ribosomal structure and biogenesis" = "#E07A5F",
  "Function unknown" = "#F2CC8F",
  "Energy production and conversion" = "#5A4B8A",
  "Nucleotide transport and metabolism" = "#7B6BB5",
  "Coenzyme transport and metabolism" = "#4F79B2",
  "Replication, recombination and repair" = "#3B5B92",
  "Transcription" = "#81B29A",
  "Amino acid transport and metabolism" = "#B7D874",
  "Carbohydrate transport and metabolism" = "#7FAE54",
  "Lipid transport and metabolism" = "#C9D86B",
  "Posttranslational modification, protein turnover, chaperones" = "#E8E287",
  "Intracellular trafficking, secretion, and vesicular transport" = "#A7B5A0",
  "Cell wall/membrane/envelope biogenesis" = "#95D5C7",
  "Signal transduction mechanisms" = "#A8DADC",
  "Secondary metabolites biosynthesis, transport and catabolism" = "#6FA8C7",
  "Cell cycle control, cell division, chromosome partitioning" = "#8FA2B7",
  "Defense mechanisms" = "#B5835A",
  "Inorganic ion transport and metabolism" = "#D96C6C",
  "Cell motility" = "#7F7F7F"
)

df_norm <- df_norm %>%
  mutate(
    Sample = factor(Sample, levels = unique(Sample)),
    COG_full = factor(cog_name, levels = cog_order)
  )

# 绘图：y使用frac（0‑1），y轴0‑100%，所有柱子严格平齐顶到100%
p_stack <- ggplot(df_norm, aes(x = Sample, y = frac, fill = COG_full)) +
  geom_col(position = position_stack(), width = 0.75) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 100),
    limits = c(0,1),
    expand = c(0,0)
  ) +
  scale_fill_manual(values = cog_colors) +
  labs(x = "Phylum", y = "Proportion (%)", fill = "COG Functional Category") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size =5),
    panel.grid = element_blank(),
    legend.text = element_text(size=5),
    legend.title = element_text(size=5)
  )

p_stack

# -------------------图2：count热图-------------------
p_heat <- ggplot(df_norm, aes(x = Sample, y = COG_full, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(low = "yellow", high = "purple") +   # 紫低黄高
  labs(x = "Phylum",
       y = "COG Functional Category",
       fill = "Relative abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank()
  )

p_heat
# 拼图
library(patchwork)
p_stack + p_heat
ggsave("result/Fig4b.tiff",p_stack + p_heat, width = 15, height = 6, dpi = 300)
ggsave("result/Fig4b.jpg",p_stack + p_heat, width = 15, height = 6, dpi = 300)
ggsave("result/Fig4b.pdf",p_stack + p_heat, width = 15, height = 6, dpi = 300)
# 首次运行先安装这些包
install.packages(c("tidyverse", "patchwork", "ggtree", "ggtreeExtra", "ape", "scales"))

# 加载包
library(tidyverse)   # 数据处理+绘图
library(patchwork)   # 子图拼接
library(ggtree)      # 进化树绘制
library(ggtreeExtra) # 进化树外圈注释
library(ape)         # 读取进化树文件
library(scales)      # 百分比格式化

library(ggplot2)
library(dplyr)
library(scales)

# ========== 1. 设置保留的防御系统类型，其余全部合并为 Others ==========
keep_defense <- c(
  "RM","Ceres","Cas","SoFIC","AbiE","MazEF","Wadjet","McrBC","Aristaios","Prometheus","CBASS","AbiAlpha"
)

df_raw <- read.csv("data/B.Phylum_defense.csv") %>%
  mutate(defense_system = ifelse(defense_system %in% keep_defense, defense_system, "Others"))
df_raw <- read.csv("data/B.Phylum_defense.csv") %>%
  mutate(
    defense_system = case_when(
      # 所有RM_Type_* 全部统一改成 RM
      grepl("^RM_Type_", defense_system) ~ "RM",
      defense_system %in% keep_defense ~ defense_system,
      TRUE ~ "Others" # Lanthiphage, SNIPE, BREX, Avs, gcu233, Other_Type_IV, DarTG, Armada 全部进Others
    )
  )

# ========== 2. 统计计数：phylum × defense_system ==========
df_count <- df_raw %>%
  count(phylum, defense_system, name = "n")

# --------------------------
# 图1：y=defense_system，fill=phylum（参考图2）
# Y轴标签：DefenseSystem(总数)
# --------------------------
top_phylum_fig1 <- df_count %>%
  group_by(phylum) %>%
  summarise(total_n = sum(n), .groups="drop") %>%
  arrange(desc(total_n)) %>%
  slice_head(n=11) %>%
  pull(phylum)

df_fig1 <- df_count %>%
  mutate(phylum = ifelse(phylum %in% top_phylum_fig1, phylum, "Others")) %>%
  group_by(defense_system) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# 计算每个defense_system总数量，拼接标签
defense_label_df <- df_fig1 %>%
  group_by(defense_system) %>%
  summarise(total = sum(n), .groups="drop") %>%
  mutate(defense_label = paste0(defense_system, "(", total, ")"))

# 标签映射 + 因子顺序
defense_level_fig1 <- c("Others","RM","Ceres","Cas","SoFIC","AbiE","MazEF","Wadjet","McrBC","Aristaios","Prometheus","CBASS","AbiAlpha")
defense_label_map <- setNames(defense_label_df$defense_label, defense_label_df$defense_system)


phylum_col <- c(
  "Actinomycetota" = "#b8d6be",
  "Acidobacteriota" = "#a4918e",
  "Bacillota" = "#e8b4b8",
  "Bacteroidota" = "#7d8ca3",
  "Desulfobacterota" = "#f2d4b7",
  "Gemmatimonadota" = "#92b6d5",
  "Planctomycetota" = "#d4a5a5",
  "Pseudomonadota" = "#c9c0bb",
  "Chloroflexota" = "#e0e6c8",
  "Cyanobacteriota" = "#917260",
  "Myxococcota" = "#b08968",
  "Others" = "#704444"
)
phylum_levels <- names(phylum_col)
df_fig1$phylum <- factor(df_fig1$phylum, levels = phylum_levels)
df_fig1$defense_system <- factor(df_fig1$defense_system, levels = defense_level_fig1)

p1 <- ggplot(df_fig1, aes(y = defense_system, x = percent, fill = phylum)) +
  geom_col(position = position_stack(), width= 0.75, color=NA, linewidth=0)+
  scale_x_continuous(labels = percent_format(), expand = c(0,0))+
  scale_fill_manual(values = phylum_col, drop = FALSE)+
  scale_y_discrete(labels = defense_label_map)+ # 替换y标签为带括号计数
  labs(x = "Number of Phylum (%)", y="Defense system", fill="Phylum")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=9))
p1
ggsave("result/Fig5a.tiff",p1, width = 6, height = 6, dpi = 300)
ggsave("result/Fig5a.jpg",p1, width = 6, height = 6, dpi = 300)
ggsave("result/Fig5a.pdf",p1, width = 6, height = 6, dpi = 300)
# --------------------------
# 图2：y=phylum，fill=defense_system（参考图3）
# Y轴标签：Phylum(总数)
# --------------------------
top_phylum_fig2 <- df_count %>%
  group_by(phylum) %>%
  summarise(total_n = sum(n), .groups="drop") %>%
  arrange(desc(total_n)) %>%
  slice_head(n=12) %>%
  pull(phylum)

df_fig2 <- df_count %>%
  mutate(phylum = ifelse(phylum %in% top_phylum_fig2, phylum, "Others")) %>%
  group_by(phylum) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# 计算每个phylum总数量，拼接标签
phylum_label_df <- df_fig2 %>%
  group_by(phylum) %>%
  summarise(total = sum(n), .groups="drop") %>%
  mutate(phylum_label = paste0(phylum, "(", total, ")"))

defense_col <- c(
  "AbiAlpha"="#e8d5c4",
  "AbiE"="#cccccc",
  "Aristaios"="#739272",
  "CBASS"="#8d9f87",
  "Cas"="#587559",
  "Ceres"="#e3c8a6",
  "MazEF"="#607890",
  "McrBC"="#b8b08d",
  "Prometheus"="#c2b4a4",
  "RM"="#b98e7c",
  "SoFIC"="#9aa8a1",
  "Wadjet"="#d6cccc",
  "Others"="#6a8596"
)
defense_level_fig2 <- names(defense_col)
df_fig2$defense_system <- factor(df_fig2$defense_system, levels = defense_level_fig2)

# Y轴phylum排序：按总数量从小到大
phylum_order_fig2 <- phylum_label_df %>%
  arrange(total) %>%
  pull(phylum)
df_fig2$phylum <- factor(df_fig2$phylum, levels = phylum_order_fig2)
phylum_label_map <- setNames(phylum_label_df$phylum_label, phylum_label_df$phylum)

p2 <- ggplot(df_fig2, aes(y = phylum, x = percent, fill = defense_system)) +
  geom_col(position = position_stack(), color=NA, linewidth=0)+
  scale_x_continuous(labels = percent_format(), expand = c(0,0))+
  scale_fill_manual(values = defense_col, drop = FALSE)+
  scale_y_discrete(labels = phylum_label_map)+ # phylum名称带括号计数
  labs(x = "Number of Defense system (%)", y="Phylum", fill="Defense system")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=9))
p2
ggsave("result/Fig5b.tiff",p2, width = 6, height = 6, dpi = 300)
ggsave("result/Fig5b.jpg",p2, width = 6, height = 6, dpi = 300)
ggsave("result/Fig5b.pdf",p2, width = 6, height = 6, dpi = 300)

library(ggplot2)
library(dplyr)

# 读入数据，假设你的数据框叫 df，两列：ID, count
df <- read.csv("data/B.defense.richness.csv")

# 1. 分组：按defense系统数量分bin
df_bin <- df %>%
  mutate(
    group = case_when(
      count == 0 ~ "0",
      count >=1 & count <=5 ~ "1~5",
      count >=6 & count <=10 ~ "6~10",
      count >=11 & count <=15 ~ "11~15",
      count >=16 & count <=20 ~ "16~20",
      count >=21 & count <=30 ~ "21~30",
      count >=31 & count <=60 ~ "31~60"
    )
  ) %>%
  count(group, name = "Genome_number")

# 设定分组顺序（保证x轴顺序不乱）
bin_order <- c("0","1~5","6~10","11~15","16~20","21~30","31~60")
df_bin$group <- factor(df_bin$group, levels = bin_order)

# 2. 绘图，复刻你示例的柱状图风格
p3 <- ggplot(df_bin, aes(x = group, y = Genome_number)) +
  geom_col(fill = "#6b7b94", width = 0.7) +
  geom_text(aes(label = Genome_number), vjust = -0.3, size = 3.8) + #柱子上方标数字
  labs(x = "Defense systems", y = "Genome number") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

p3
ggsave("result/Fig5c.tiff",p3, width = 4, height = 6, dpi = 300)
ggsave("result/Fig5c.jpg",p3, width = 4, height = 6, dpi = 300)
ggsave("result/Fig5c.pdf",p3, width = 4, height = 6, dpi = 300)
library(ggplot2)
library(dplyr)
library(patchwork)
provirus_df <- read.csv("data/checkv.1750.csv")
# ====================== 数据统计 ======================
# 1. 水平条形图：phylum计数
phylum_count <- provirus_df %>%
  count(phylum, name="count") %>%
  arrange(count) %>% # 水平条形图y轴从小到大
  mutate(phylum = factor(phylum, levels=phylum))

# 配色
phylum_col <- c(
  "Pseudomonadota"="#ffdd70",
  "Bacteroidota"="#c8e8bc",
  "Actinomycetota"="#b078bc",
  "Bacillota"="#f8c6e0",
  "Gemmatimonadota"="#b8dd70",
  "Chloroflexota"="#ffaa55",
  "Planctomycetota"="#77aacc",
  "Myxococcota"="#f2705c",
  "Bacillota_I"="#abcdef",
  "Acidobacteriota"="#fff299",
  "Bacillota_A"="#e3a6dd",
  "Halobacteriota"="#78c8bc",
  "Thermoplasmatota"="#e6cccc",
  "Others"="#d6d6d6"
)

comp_col <- c(
  "Complete"="#79b873",
  "High-quality"="#b28868",
  "Medium-quality"="#e6cc80"
)

# ====================== 绘图 ======================
# 左：水平条形图
p_bar <- ggplot(phylum_count, aes(y=phylum, x=count, fill=phylum)) +
  geom_col(width=0.7, color=NA)+
  geom_text(aes(label=count), hjust=-0.2, size=3.5)+
  scale_fill_manual(values=phylum_col, guide="none")+
  labs(x="Count", y="Phylum")+
  theme_bw()+
  theme(panel.grid = element_blank())
p_bar
ggsave("result/Fig5d.tiff",p_bar, width = 6, height = 6, dpi = 300)
ggsave("result/Fig5d.jpg",p_bar, width = 6, height = 6, dpi = 300)
ggsave("result/Fig5d.pdf",p_bar, width = 6, height = 6, dpi = 300)
library(ggplot2)
library(dplyr)
library(patchwork)

# ====================== 1. 统计数据 + 合并低丰度为Others ======================
# 1.1 phylum计数，合并低丰度门类
phylum_count <- provirus_df %>%
  count(phylum, name="count") %>%
  # 保留Top10，其余全部合并成Others
  mutate(phylum = ifelse(rank(-count) <= 15, phylum, "Others")) %>%
  group_by(phylum) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  arrange(count) %>%
  mutate(phylum = factor(phylum, levels = phylum))

# 1.2 饼图I 
checkv_count <- provirus_df %>%
  count(checkv_quality, name="cnt") %>%
  mutate(percent = round(cnt / sum(cnt)*100,2))

# 1.3 饼图II host type：virus-bacteria(301), virus-archaea(4)
host_type_df <- data.frame(
  host = c("virus–bacteria", "Others"),
  cnt = c(190, 1560)
) %>% mutate(percent = round(cnt/sum(cnt)*100,2))

# 1.4 饼图III category
cat_df <- data.frame(
  group = c("virus-archaea", "Others"),
  cnt = c(4,174)
) %>% mutate(percent = round(cnt/sum(cnt)*100,2))

# ====================== 配色（新增Others灰色） ======================
phylum_col <- c(
  "Pseudomonadota"="#ffdd70",
  "Bacteroidota"="#c8e8bc",
  "Actinomycetota"="#b078bc",
  "Bacillota"="#f8c6e0",
  "Gemmatimonadota"="#b8dd70",
  "Chloroflexota"="#ffaa55",
  "Planctomycetota"="#77aacc",
  "Myxococcota"="#f2705c",
  "Bacillota_I"="#abcdef",
  "Acidobacteriota"="#fff299",
  "Bacillota_A"="#e3a6dd",
  "Halobacteriota"="#78c8bc",
  "Bacillota_D"="#e6cccc",
  "Others"="#d6d6d6"
)
checkv_col <- c(
  "Complete"="#79b873",
  "High-quality"="#b28868",
  "Medium-quality"="#e6cc80"
)
host_col <- c("virus–bacteria"="#82a8d8", "Others"="#e6edf7")
cat_col <- c("virus-archaea"="#b4b0d8", "Others"="#eeedf5")

# ====================== 绘图 ======================
# 左：水平条形图
p_bar <- ggplot(phylum_count, aes(y=phylum, x=count, fill=phylum)) +
  geom_col(width=0.7, color=NA)+
  geom_text(aes(label=count), hjust=-0.2, size=3.5)+
  scale_fill_manual(values=phylum_col, guide="none")+
  labs(x="Count", y="Phylum")+
  theme_bw()+
  theme(panel.grid = element_blank())
p_bar
# 饼图 I completeness
p_pie1 <- ggplot(checkv_count, aes(x="", y=cnt, fill=checkv_quality)) +
  geom_col()+
  coord_polar("y", start=0)+
  geom_text(aes(label=paste0(percent,"%")), position=position_stack(vjust=0.5), size=3.2)+
  scale_fill_manual(values=checkv_col)+
  labs(title="I")+
  theme_void()
p_pie1
# 饼图 II host type
p_pie2 <- ggplot(host_type_df, aes(x="", y=cnt, fill=host)) +
  geom_col()+
  coord_polar("y", start=0)+
  geom_text(aes(label=paste0(percent,"%")), position=position_stack(vjust=0.5), size=3.2)+
  scale_fill_manual(values=host_col)+
  labs(title="II")+
  theme_void()
p_pie2
# 饼图 III category
p_pie3 <- ggplot(cat_df, aes(x="", y=cnt, fill=group)) +
  geom_col()+
  coord_polar("y", start=0)+
  geom_text(aes(label=paste0(percent,"%")), position=position_stack(vjust=0.5), size=3.2)+
  scale_fill_manual(values=cat_col)+
  labs(title="III")+
  theme_void()
p_pie3
# ====================== 拼图布局 ======================
right_panel <- (p_pie1 / p_pie2 / p_pie3)
right_panel
ggsave("result/Fig5dd.tiff",right_panel, width = 3, height = 5, dpi = 300)
ggsave("result/Fig5dd.jpg",right_panel, width = 3, height = 5, dpi = 300)
ggsave("result/Fig5dd.pdf",right_panel, width = 3, height = 5, dpi = 300)

## ================= e图 环形进化树+外圈注释 =================
tree <- read.tree("data/gtdbtk_iqtree2_result.treefile") #替换成你的nwk树文件
tree <- ladderize(tree)
library(ggtreeExtra)
library(ggplot2)
library(ggtree)
library(ape)
group_color <- c(
  "Acidobacteriota"="#b6d7a8",
  "Actinomycetota"="#b4a79b",
  "Bacillota"="#f4cccc",
  "Bacteroidota"="#6fa8dc",
  "Chloroflexota"="#e5efd9",
  "Cyanobacteriota"="#c27ba0",
  "Desulfobacterota"="#f8cbad",
  "Gemmatimonadota"="#76a5af",
  "Myxococcota"="#e6ffff",
  "Others"="#dddddd",
  "Planctomycetota"="#ea9999",
  "Pseudomonadota"="#d0c9c0",
  "Verrucomicrobiota"="#996644"
)
tree.plot <- ggtree(tree, layout = "fan", right = TRUE, open.angle = 15, size = 0.1) %<+% metadata +
  geom_aline(aes(color = Phylum..gtdb.), linetype = "solid", size = 0.2, show.legend = TRUE) +
  scale_color_manual(values = group_color, name = "Phylum") +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
    )
  ) +
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.5, "mm")
  )
tree.plot
# 按id分组，统计每个id出现多少次
meta_count <- provirus_df %>%
  group_by(id) %>%
  summarise(count = n(), .groups = "drop")

# 查看结果，11-A.bin.5 的count就等于2
meta_count
library(dplyr)
# 查看count的最大、最小
meta_count %>% summarise(
  min_count = min(count),
  max_count = max(count)
)

# 查看哪一行是最大值、最小值
meta_count %>% slice_max(count, n=1)
meta_count %>% slice_min(count, n=1)
library(ggnewscale)

tree.plot <- ggtree(tree, layout = "fan", right = TRUE, open.angle = 15, size = 0.1) %<+% metadata +
  # ========= 先放geom_fruit！！外圈柱状图放最前面 =========
new_scale_fill() +
  geom_fruit(
    data = meta_count,
    geom = geom_col,
    mapping = aes(y = id, x = count),
    pwidth = 0.2,
    orientation = "y",
    fill = "#222222"
  ) +
  # ========= 再画树的aline线和其他图层 =========
geom_aline(aes(color = Phylum..gtdb.), linetype = "solid", size = 0.2, show.legend = TRUE) +
  scale_color_manual(values = group_color, name = "Phylum") +
  guides(
    colour = guide_legend(
      title.position = "top",
      title.theme = element_text(face = "bold", size = 12, margin = margin(b = 20))
    )
  ) +
  theme(
    legend.position = "left",
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.spacing.y = unit(1.5, "mm")
  )
tree.plot
ggsave("result/Fig5e.tiff",tree.plot, width = 6, height = 5, dpi = 300)
ggsave("result/Fig5e.jpg",tree.plot, width = 6, height = 5, dpi = 300)
ggsave("result/Fig5e.pdf",tree.plot, width = 6, height = 5, dpi = 300)

library(ggplot2)
library(dplyr)
library(scales)

# ========== 1. 设置保留的防御系统类型，其余全部合并为 Others ==========
keep_defense <- c(
  "Ceres","RM","SoFIC","AbiE","Prometheus","Cas","AbiAlpha","AbiU", "ScoMcrA", "DpnI", "Aristaios","pAgo"
)


df_raw <- read.csv("data/A.Phylum_defense.csv") %>%
  mutate(
    defense_system = case_when(
      # 所有RM_Type_* 全部统一改成 RM
      grepl("^RM_Type_", defense_system) ~ "RM",
      defense_system %in% keep_defense ~ defense_system,
      TRUE ~ "Others" # Lanthiphage, SNIPE, BREX, Avs, gcu233, Other_Type_IV, DarTG, Armada 全部进Others
    )
  )

# ========== 2. 统计计数：phylum × defense_system ==========
df_count <- df_raw %>%
  count(phylum, defense_system, name = "n")

# --------------------------
# 图1：y=defense_system，fill=phylum（参考图2）
# Y轴标签：DefenseSystem(总数)
# --------------------------
top_phylum_fig1 <- df_count %>%
  group_by(phylum) %>%
  summarise(total_n = sum(n), .groups="drop") %>%
  arrange(desc(total_n)) %>%
  slice_head(n=11) %>%
  pull(phylum)

df_fig1 <- df_count %>%
  mutate(phylum = ifelse(phylum %in% top_phylum_fig1, phylum, "Others")) %>%
  group_by(defense_system) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# 计算每个defense_system总数量，拼接标签
defense_label_df <- df_fig1 %>%
  group_by(defense_system) %>%
  summarise(total = sum(n), .groups="drop") %>%
  mutate(defense_label = paste0(defense_system, "(", total, ")"))

# 标签映射 + 因子顺序
defense_level_fig1 <- c("Ceres","Others","RM","SoFIC","AbiE","Prometheus","Cas","AbiAlpha","AbiU", "ScoMcrA", "DpnI", "Aristaios","pAgo")
defense_label_map <- setNames(defense_label_df$defense_label, defense_label_df$defense_system)
phylum_col <- c(
  "Halobacteriota" = "#4a5f76",
  "Thermoplasmatota" = "#4c6b57",
  "Thermoproteota" = "#d9b38c",
  "Nanoarchaeota"="#b0cfc7",
  "Nanohaloarchaeota"="#aaabc2"
)

phylum_levels <- names(phylum_col)
df_fig1$phylum <- factor(df_fig1$phylum, levels = phylum_levels)
df_fig1$defense_system <- factor(df_fig1$defense_system, levels = defense_level_fig1)

p1 <- ggplot(df_fig1, aes(y = defense_system, x = percent, fill = phylum)) +
  geom_col(position = position_stack(), width= 0.75, color=NA, linewidth=0)+
  scale_x_continuous(labels = percent_format(), expand = c(0,0))+
  scale_fill_manual(values = phylum_col, drop = FALSE)+
  scale_y_discrete(labels = defense_label_map)+ # 替换y标签为带括号计数
  labs(x = "Number of Phylum (%)", y="Defense system", fill="Phylum")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=9))
p1
ggsave("result/sFig11a.tiff",p1, width = 6, height = 6, dpi = 300)
ggsave("result/sFig11a.jpg",p1, width = 6, height = 6, dpi = 300)
ggsave("result/sFig11a.pdf",p1, width = 6, height = 6, dpi = 300)

top_phylum_fig2 <- df_count %>%
  group_by(phylum) %>%
  summarise(total_n = sum(n), .groups="drop") %>%
  arrange(desc(total_n)) %>%
  slice_head(n=12) %>%
  pull(phylum)

df_fig2 <- df_count %>%
  mutate(phylum = ifelse(phylum %in% top_phylum_fig2, phylum, "Others")) %>%
  group_by(phylum) %>%
  mutate(percent = n / sum(n)) %>%
  ungroup()

# 计算每个phylum总数量，拼接标签
phylum_label_df <- df_fig2 %>%
  group_by(phylum) %>%
  summarise(total = sum(n), .groups="drop") %>%
  mutate(phylum_label = paste0(phylum, "(", total, ")"))

defense_col <- c(
  "pAgo"="#e8d5c4",
  "Aristaios"="#cccccc",
  "DpnI"="#739272",
  "ScoMcrA"="#8d9f87",
  "AbiU"="#587559",
  "AbiAlpha"="#e3c8a6",
  "Cas"="#607890",
  "Prometheus"="#b8b08d",
  "AbiE"="#c2b4a4",
  "SoFIC"="#b98e7c",
  "RM"="#9aa8a1",
  "Ceres"="#6a8596",
  "Others"="#d0c9c0"
)
defense_level_fig2 <- names(defense_col)
df_fig2$defense_system <- factor(df_fig2$defense_system, levels = defense_level_fig2)

# Y轴phylum排序：按总数量从小到大
phylum_order_fig2 <- phylum_label_df %>%
  arrange(total) %>%
  pull(phylum)
df_fig2$phylum <- factor(df_fig2$phylum, levels = phylum_order_fig2)
phylum_label_map <- setNames(phylum_label_df$phylum_label, phylum_label_df$phylum)

p2 <- ggplot(df_fig2, aes(y = phylum, x = percent, fill = defense_system)) +
  geom_col(position = position_stack(), color=NA, linewidth=0)+
  scale_x_continuous(labels = percent_format(), expand = c(0,0))+
  scale_fill_manual(values = defense_col, drop = FALSE)+
  scale_y_discrete(labels = phylum_label_map)+ # phylum名称带括号计数
  labs(x = "Number of Defense system (%)", y="Phylum", fill="Defense system")+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=9))
p2
ggsave("result/sFig11b.tiff",p2, width = 6, height = 6, dpi = 300)
ggsave("result/sFig11b.jpg",p2, width = 6, height = 6, dpi = 300)
ggsave("result/sFig11b.pdf",p2, width = 6, height = 6, dpi = 300)
library(ggplot2)
library(dplyr)

# 读入数据，假设你的数据框叫 df，两列：ID, count
df <- read.csv("data/A.defense.richness.csv")

# 1. 分组：按defense系统数量分bin
df_bin <- df %>%
  mutate(
    group = case_when(
      count == 0 ~ "0",
      count >=1 & count <=5 ~ "1~5",
      count >=6 & count <=10 ~ "6~10",
      count >=11 & count <=15 ~ "11~15",
      count >=16 & count <=20 ~ "16~20",
      count >=21 & count <=30 ~ "21~30",
      count >=31 & count <=60 ~ "31~60"
    )
  ) %>%
  count(group, name = "Genome_number")

# 设定分组顺序（保证x轴顺序不乱）
bin_order <- c("0","1~5","6~10","11~15","16~20","21~30","31~60")
df_bin$group <- factor(df_bin$group, levels = bin_order)

# 2. 绘图，复刻你示例的柱状图风格
p3 <- ggplot(df_bin, aes(x = group, y = Genome_number)) +
  geom_col(fill = "#aaabc2", width = 0.7) +
  geom_text(aes(label = Genome_number), vjust = -0.3, size = 3.8) + #柱子上方标数字
  labs(x = "Defense systems", y = "Genome number") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

p3
ggsave("result/sFig11c.tiff",p3, width = 4, height = 6, dpi = 300)
ggsave("result/sFig11c.jpg",p3, width = 4, height = 6, dpi = 300)
ggsave("result/sFig11c.pdf",p3, width = 4, height = 6, dpi = 300)
library(patchwork)

final_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(
    tag_levels = "a", # 自动A,B,C
    tag_prefix = "", tag_suffix = ""
  ) & theme(plot.tag.position = c(0.05, 0.95)) #标签左上角
final_plot
ggsave("result/sFig11.tiff",final_plot, width = 18, height = 6, dpi = 300)
ggsave("result/sFig11.jpg",final_plot, width = 18, height = 6, dpi = 300)
ggsave("result/sFig11.pdf",final_plot, width = 18, height = 6, dpi = 300)

