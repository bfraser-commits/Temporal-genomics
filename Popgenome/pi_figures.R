library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

pi <- data.frame(fread("Documents/FIBR_temporal/popgenome/ALL.50kb.pi.popgenome.out"))
head(pi)
colnames(pi) <- c("chrom","window","window_start","window_end","GH","CA13",
                  "CA18","TA13","TA18","LL13","LL18","UL13","UL18")


long_pi <- pi %>%
  pivot_longer(
    cols = starts_with(c("CA", "TA", "LL", "UL")),
    names_to = "pop_year",
    values_to = "PI"
  ) %>%
  mutate(
    # Extract population
    pop = case_when(
     # pop_year == "GH" ~ "GH",
      TRUE ~ str_extract(pop_year, "[A-Z]+")
    ),
    # Extract year
    year = case_when(
     # pop_year == "GH" ~ "GH",
      str_detect(pop_year, "13") ~ "13",
      str_detect(pop_year, "18") ~ "18",
      TRUE ~ NA_character_
    ),
    pop = factor(pop, levels = c("UL", "LL", "CA", "TA"))
  ) %>%
  filter(!is.na(year))

head(long_pi)
table(long_pi$pop)

GH_plot <- ggplot(pi, aes(x= "GH", y=GH))+
  geom_violin()+
  geom_boxplot(
  width = 0.1,
  outlier.shape = NA,
  alpha = 0.5,
  position = position_dodge(width = 0.8)
)   +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  ) +
  xlab("") +
  ylab("PI")

table(long_pi$year)


POP_pi_plot <- ggplot(long_pi, aes(x = year, y = PI, fill = pop)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7,
    position = position_dodge(width = 0.9)
  ) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA,
    alpha = 0.5,
    position = position_dodge(width = 0.9)
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = c(
    UL = "#CD0000",
    LL = "#5D478B",
    CA = "#FFB90F",
    TA = "#FF7F00"
  )) +
  xlab("Year") +
  ylab("PI")+
  coord_cartesian(ylim = c(0, 0.012))

plot_grid(GH_plot, POP_pi_plot, ncol = 2, align = "h")


GH_data <- pi %>%
  select(GH) %>%
  mutate(pop = "GH")

# Now plot GH as a "population" like the others
GH_plot <- ggplot(GH_data, aes(x = pop, y = GH, fill = pop)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.9) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"  # Remove legend for single group
  ) +
  scale_fill_manual(values = c(GH = "#00008B")) +
  xlab("") +
  ylab("PI")+
  coord_cartesian(ylim = c(0, 0.012))

# Combine plots with appropriate relative widths
violin_plot <- plot_grid(GH_plot, POP_pi_plot, ncol = 2, align = "h", rel_widths = c(1, 3))

#range(pi$GH, long_pi$PI, na.rm = TRUE)

png("Documents/FIBR_temporal/popgenome/pi_violin_plot.png", width=8, height=4, units = "in", res = 300)
violin_plot
dev.off()

##mean figures
head(pi)
#GH_13_cols <- c("GH_CA13", "GH_TA13", "GH_LL13", "GH_UL13")
pi_subset <- pi[, grepl("^(GH|CA|TA|LL|UL)", names(pi))]
head(pi_subset)

pi_means <- colMeans(pi_subset, na.rm = TRUE)
pi_se <- sapply(pi_subset, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
pi_means

pi_pop <- data.frame("pop" = c("GH","CA","CA", "TA","TA","LL", "LL", "UL","UL"))
pi_pop$mean <- as.numeric(pi_means)
pi_pop$SE <- as.numeric(pi_se)
pi_pop$year <- c("GH", "2013", "2018", "2013", "2018", "2013", "2018", "2013", "2018")

mean_pi_plot <- ggplot(pi_pop, aes(x = year, y = mean, colour = pop)) +
  #geom_point(size = 2)+
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  # geom_errorbar(
  #   aes(ymin = mean - SE, ymax = mean + SE),
  #   width = 0.2,
  #   position = position_dodge(width = 0.2),
  # )+
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.position = "none",
    plot.margin = margin(c(0.5, 0.7, 0.1, 0.1), unit = "cm")
  ) +
  scale_colour_manual(values = c(GH = "#00008B", CA = "#FFB90F", TA = "#FF7F00", LL = "#5D478B", UL = "#CD0000")) +
  xlab(expression("Year")) +
  ylab(expression("PI")) 
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 0.06), breaks = c(0, 0.02, 0.04, 0.06),
   #                  labels = c(0, 0.02, 0.04, 0.06))

mean_pi_plot
