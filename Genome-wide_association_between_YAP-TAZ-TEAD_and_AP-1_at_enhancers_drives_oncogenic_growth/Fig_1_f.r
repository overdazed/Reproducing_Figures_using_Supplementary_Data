library(rtracklayer)
library(here)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(ggplot2)
library(scales)

TAZ_peaks <- import(here("data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peak<- import(here("data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))


# -----------------------------
# 1. Get TSS
# -----------------------------
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
tss_gr <- promoters(hg38_transcripts, upstream = 0, downstream = 1)

# -----------------------------
# 2. Compute distances for YAP, TAZ, TEAD4
# -----------------------------
YAP_dist <- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist <- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist <- mcols(distanceToNearest(TEAD4_peak, tss_gr))$distance

# Put them in a single data frame
tss_distance_df <- bind_rows(
  data.frame(factor = "YAP", distance = YAP_dist),
  data.frame(factor = "TAZ", distance = TAZ_dist),
  data.frame(factor = "TEAD4", distance = TEAD4_dist)
)

# -----------------------------
# 3. Compute distances for overlapping peaks
# -----------------------------

# Create a GRanges object for the overlap of all three
YAP_TAZ_overlap <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
YAP_TAZ_TEAD_overlap <- subsetByOverlaps(YAP_TAZ_overlap, TEAD4_peak)

Overlap_dist <- mcols(distanceToNearest(YAP_TAZ_TEAD_overlap, tss_gr))$distance
overlap_df <- data.frame(factor = "Overlap", distance = Overlap_dist)

# Combine with existing data frame
tss_distance_all <- bind_rows(tss_distance_df, overlap_df)

# length(YAP_TAZ_TEAD_overlap)  # number of overlapping peaks
# length(Overlap_dist)           # should match the number of rows in overlap_df


# -----------------------------
# 4. Categorize distances
# -----------------------------
tss_distance_all <- tss_distance_all %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >= 1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <= 100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  ))

# -----------------------------
# 5. Count per category and total
# -----------------------------
counts_all <- tss_distance_all %>%
  group_by(factor, category) %>%
  count()

totals_all <- tss_distance_all %>%
  group_by(factor) %>%
  count(name = "total")

# Merge and compute percentage
merged_all <- left_join(counts_all, totals_all, by = "factor") %>%
  mutate(Percentage = n / total * 100)

# -----------------------------
# 6. Set factor levels for correct bar order
# -----------------------------
merged_all$factor <- factor(merged_all$factor, 
                            levels = c("YAP", "TAZ", "TEAD4", "Overlap"))
merged_all$category <- factor(merged_all$category,
                              levels = c("<1kb", "1-10kb", "10-100kb", ">100kb"))

# -----------------------------
# 7. Plot stacked bar
# -----------------------------
ggplot(merged_all, aes(x = factor, y = Percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Distance to TSS", x = "Group", y = "Percentage") +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    expand = c(0, 0)  # remove padding at bottom and top
  ) +
  scale_fill_manual(values = c("#EF3E2B", "#F16161", "#F59595", "#FAD1C8")) +
  theme_classic(base_size = 14)

# cat("Number of YAP peaks:", length(YAP_peaks), "\n")
# cat("Number of TAZ peaks:", length(TAZ_peaks), "\n")
# cat("Number of TEAD4 peaks:", length(TEAD4_peak), "\n")
# cat("Number of triple overlaps:", length(YAP_TAZ_TEAD_overlap), "\n")

#ggsave("fig_1_f_overlap.png", width = 10, height = 10, units = "in")