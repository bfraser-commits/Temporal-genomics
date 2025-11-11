library(tidyr)
library(dplyr)
library(ggplot2)

FST <- data.frame(fread("Documents/FIBR_temporal/popgenome/ALL.50kb.fst.popgenome.out"))
head(FST)
colnames(FST) <- c("chrom", "window", "window_start","window_end","GH_CA13", "GH_CA18", 
                   "GH_TA13", "GH_TA18", "GH_LL13", "GH_LL18", "GH_UL13", "GH_UL18", 
                   "CA13_CA18","CA13_TA13","CA13_TA18","CA13_LL13","CA13_LL18","CA13_UL13",
                   "CA13_UL18","CA18_TA13","CA18_TA18","CA18_LL13","CA18_LL18","CA18_UL13",
                   "CA18_UL18","TA13_TA18","TA13_LL13","TA13_LL18","TA13_UL13","TA13_UL18",
                   "TA18_LL13","TA18_LL18","TA18_UL13","TA18_UL18","LL13_LL18","LL13_UL13",
                   "LL13_UL18","LL18_UL13","LL18_UL18","UL13_UL18")

summary(FST)
GH_13_cols <- c("GH_CA13", "GH_TA13", "GH_LL13", "GH_UL13")
GH_13_subset <- FST[ , GH_13_cols]
GH_13_means <- colMeans(GH_13_subset, na.rm = TRUE)
GH_13_se <- sapply(GH_13_subset, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

head(FST)


FST_pop <- data.frame("pop" = c("CA", "TA", "LL", "UL"))
FST_pop$mean <- as.numeric(GH_13_means)
FST_pop$SE <- as.numeric(GH_13_se)
FST_pop$year <- c("2013", "2013", "2013", "2013")

GH_18_cols <- c("GH_CA18", "GH_TA18", "GH_LL18", "GH_UL18")
GH_18_subset <- FST[ , GH_18_cols]
GH_18_means <- colMeans(GH_18_subset, na.rm = TRUE)
GH_18_se <- sapply(GH_18_subset, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))


FST_pop_18 <- data.frame("pop" = c("CA", "TA", "LL", "UL"))
FST_pop_18$mean <- as.numeric(GH_18_means)
FST_pop_18$SE <- as.numeric(GH_18_se)
FST_pop_18$year <- c("2018", "2018", "2018", "2018")
FST_pop_18

Inter_pop <- c("CA13_CA18", "TA13_TA18", "LL13_LL18", "UL13_UL18")
Inter_pop_subset <- FST[ , Inter_pop]
Inter_pop_means <- colMeans(Inter_pop_subset, na.rm = TRUE)
Inter_pop_subset_se <- sapply(Inter_pop_subset, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))


Inter_pop_df <- data.frame("pop" = c("CA", "TA", "LL", "UL"))
Inter_pop_df$mean <- as.numeric(Inter_pop_means)
Inter_pop_df$SE <- as.numeric(Inter_pop_subset_se)
Inter_pop_df$year <- c("2018", "2018", "2018", "2018")
Inter_pop_df
FST_pop

Final_df <- data.frame(rbind(FST_pop, FST_pop_18))

Final_df
Final_df$pop <- factor(Final_df$pop, levels = c("LL", "UL", "CA", "TA"))
Final_df$year <- factor(Final_df$year)

mean_FST_plot <- ggplot(Final_df, aes(x = year, y = mean, colour = pop)) +
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
  scale_colour_manual(values = c(CA = "#FFB90F", TA = "#FF7F00", LL = "#5D478B", UL = "#CD0000")) +
  xlab(expression("Year")) +
  ylab(expression("FST")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.06), breaks = c(0, 0.02, 0.04, 0.06),
                     labels = c(0, 0.02, 0.04, 0.06))

mean_FST_plot
pdf("Documents/FIBR_temporal/popgenome/FST_mean_plot.pdf", width=5, height=4)
mean_FST_plot
dev.off()

##try violin plots instead
#first get df that has pop column 1, year column 2, and FST column 3
raw_data <- data.frame(cbind(GH_13_subset, GH_18_subset))
head(raw_data)  

long_df <- raw_data %>%
  pivot_longer(
    cols = everything(),  # pivot all columns
    names_to = c(".value", "pop_year"), # temporarily keep GH prefix, then pop_year combined
    names_pattern = "(GH)_(.*)"  # split into 'GH' and the rest (pop_year)
  ) %>%
  # Now separate 'pop_year' into pop and year parts
  separate(pop_year, into = c("pop", "year_suffix"), sep = -2) %>%
  mutate(
    year = as.integer(paste0("20", year_suffix)),   # convert 13 -> 2013, 18 -> 2018
    pop = factor(pop, levels = c("UL", "LL", "CA", "TA"))  # order populations
  ) %>%
  rename(value = GH) %>%   # rename the measurement column
  select(pop, year, value) # keep relevant columns

summary(long_df)


violin_Fst_plot <- ggplot(long_df, aes(x = factor(year), y = value, fill = pop)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c(UL = "#CD0000", LL = "#5D478B", CA = "#FFB90F", TA = "#FF7F00")) +
  scale_colour_manual(values = c(UL = "#CD0000", LL = "#5D478B", CA = "#FFB90F", TA = "#FF7F00")) +
  xlab("Year") +
  ylab("FST")

png("Documents/FIBR_temporal/popgenome/FST_violin_plot.png", width = 8, height=4, units = "in", res = 300)
violin_Fst_plot
dev.off()

dodge_width <- 0.9  # You can adjust this to your liking

###with boxplots###
violin_Fst_plot <- ggplot(long_df, aes(x = factor(year), y = value, fill = pop)) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7,
    position = position_dodge(width = dodge_width),
    width = 0.9
  ) +
  geom_boxplot(
    width = 0.1,
    alpha = 0.5,
    outlier.shape = NA,
    position = position_dodge(width = dodge_width)
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.position = "right"
  ) +
  scale_fill_manual(values = c(
    UL = "#CD0000", LL = "#5D478B",
    CA = "#FFB90F", TA = "#FF7F00"
  )) +
  xlab("Year") +
  ylab("FST")

png("Documents/FIBR_temporal/popgenome/FST_violin_plot_box.png", width = 8, height=4, units = "in", res = 300)
violin_Fst_plot
dev.off()
