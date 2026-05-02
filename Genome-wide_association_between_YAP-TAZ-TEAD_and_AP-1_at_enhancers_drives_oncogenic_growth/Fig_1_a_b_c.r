# Big picture: what this analysis actually is

# You are reproducing Figure 1 a, b, c from a ChIP-seq paper.

# Biological question:
  
#  - Do YAP and TAZ bind the same genomic regions?
  
#  - Do those shared regions also recruit TEAD4?
  
#  - Where is TEAD4 positioned relative to YAP/TAZ peak summits?
  
# Everything in the script is geometry on the genome.
# No statistics yet. Just intervals and distances.

# Use GenomicRanges in R
# https://pubmed.ncbi.nlm.nih.gov/26258633/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66081

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("ggplot2")

# Let’s use R instead with the `GenomicRanges` package. Read my [blog post](https://divingintogeneticsandgenomics.com/post/genomic-interval/) on why we need to learn how to use it.


library(GenomicRanges)
library(rtracklayer) # for reading in bed file, bioConduct package, use to read in the peak files
library(here) # to give full path of the file

#=======================================================================

# ----- 1. Loading peak files → GenomicRanges ----- #

# What happens here:

#  - rtracklayer::import() reads a BED-like file

#  - It creates a GRanges object

#  - Each row = one peak

#  - Each peak has:
#    - seqnames (chromosome)
#    - ranges (start, end)
#    - metadata columns (score, signalValue, etc.)

#-----------------------------------------------------------------------

# TAZ_peaks <- import(here("1_reproduce_figures_using_Supplementary_Data/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
TAZ_peaks <- import("N:/Bioinformatik/Programming/Projects/1_reproduce_figures_using_Supplementary_Data/data/fastq/TAZ_peak/TAZ_peaks.narrowPeak")

# file.exists(here("data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
# list.files(here("data/fastq/TAZ_peak"))

# list.files(
#   "N:/Bioinformatik/Programming/Projects",
#   pattern = "narrowPeak",
#   recursive = TRUE,
#   full.names = TRUE
# )

YAP_peaks <- import("N:/Bioinformatik/Programming/Projects/1_reproduce_figures_using_Supplementary_Data/data/fastq/YAP_peak/YAP_peaks.narrowPeak")
#YAP_peaks<- import(here("1_reproduce_figures_using_Supplementary_Data/data/fastq/YAP_peak/YAP_peaks.narrowPeak"))

TAZ_peaks # gives you the GenomicRanges object, 10719
YAP_peaks # gives you the GenomicRanges object,  9807

#=======================================================================

# ----- 2. Overlap logic: YAP vs TAZ ----- #

# Important conceptual point:
#  - subsetByOverlaps(A, B)
#    → “keep A peaks that overlap any peak in B”

# That means:
#  - length(A ∩ B) is not symmetric
#  - One peak can overlap multiple peaks in the other set

# That’s why you get:
#  - 7154 TAZ peaks overlapping YAP
#  - 7164 YAP peaks overlapping TAZ


# Hidden assumption:
#  - “Overlap” means ≥1 bp
#  - No distance threshold
#  - No summit logic yet

# ----------------------------------------------------------------------

# We have `r length(YAP_peaks)` YAP peaks and `r length(TAZ_peaks)` TAZ peaks.
# How many of them overlap?
TAZ_overlap_YAP_peaks <- subsetByOverlaps(TAZ_peaks, YAP_peaks)
length(TAZ_overlap_YAP_peaks) # 7154

TAZ_overlap_YAP_peaks # look at it

YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)
length(YAP_overlap_TAZ_peaks) # 7164

YAP_overlap_TAZ_peaks # look at it

# They are not exactly the same because some peaks overlap with multiple peaks
# from the other set.

# so we have 7154 out of 10719 TAZ peaks overlapping with YAP1 peaks.
# (Note, this is the same number we got from bedtools intersect).

# and we have 7164 out of 9807 YAP1 peaks overlapping with TAZ peaks.

# The venn-diagram needs a common number in the intersection.
# How do we deal with it? There are different decisions you can make
# and it does not affect the conclusion of the figure:
#  most of the TAZ and YAP1 peaks overlap!

# We can just use the number of YAP1 peaks that overlap with TAZ as the intersection.

#=======================================================================

# ----- 3. The Venn diagram problem -----

# There are many packages to make a venndiagram. I use [Vennerable.](https://github.com/js229/Vennerable)

# devtools::install_github("js229/Vennerable")

#install.packages("devtools")
#  library(devtools)

# install.packages("BiocManager")
#BiocManager::install(version = "3.19")
#  BiocManager::install(version = "3.22")

# BiocManager::install(c("graph", "RBGL"))

# Sys.which("make")

# devtools::install_github("js229/Vennerable")
# remotes::install_github("js229/Vennerable")
library(Vennerable)


# Venn diagrams want one number for the intersection.

# Reality:
#  - Overlaps are many-to-many
#  - Biology doesn’t care about your plotting library


n_YAP <- length(YAP_peaks)  # Total number of peaks for each Set
n_TAZ <- length(TAZ_peaks)  # Total number of peaks for each Set


# This is a reasonable simplification, and the paper explicitly allows it.
# Key idea:
#   The conclusion (“most peaks overlap”) is robust to the exact counting choice.
n_overlap <- length(YAP_overlap_TAZ_peaks)

venn_data <- Venn(SetNames = c("αYAP peaks", "αTAZ peaks"),
                  Weight = c(
                    # 1 = only in set 1, but not in set 2
                    "10" = n_YAP, # Unique to A
                    # 2 = only in set 2, but not in set 1
                    "01" = n_TAZ, # Unique to B
                    # common to both sets
                    "11" = n_overlap         # Intersection
                  ))


# Plot the Venn diagram
# plot(venn_data)

plot_fig1a <- plot(venn_data)
# You can also use makeVennDiagram in the ChIPpeakAnno package.

# In our case, we already have the number of the two sets and the intersection, 
# so we used vennerable.

# Take a look at [ggVennDiagram](https://github.com/gaospecial/ggVennDiagram) too

##### Figure 1 b
# We can easily make Figure 1b now that we have some foundations.

# TEAD4_peak <- import(here("1_reproduce_figures_using_Supplementary_Data/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))
TEAD4_peak <- import(here("N:/Bioinformatik/Programming/Projects/1_reproduce_figures_using_Supplementary_Data/data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

#=======================================================================

# ----- 4. Adding TEAD4: second-order overlap -----

# Now the logic stack is:
#  1. Start with YAP peaks
#  2. Keep only those overlapping TAZ
#  3. From those, keep only those overlapping TEAD4

# This answers:
# “Among YAP/TAZ shared sites, how many recruit TEAD4?”

YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)

n_YAP_TAZ <- length(YAP_overlap_TAZ_peaks)  # Total peaks 
n_TEAD4 <- length(TEAD4_peak)  # Total peaks 
n_overlap2 <- length(YAP_overlap_TAZ_peaks_overlap_TEAD4)

venn_data2 <- Venn(SetNames = c("αYAP/TAZ peaks", "αTEAD4 peaks"),
                  Weight = c(
                    "10" = n_YAP_TAZ, # Unique to A
                    "01" = n_TEAD4, # Unique to B
                    "11" = n_overlap2        # Intersection
                  ))

# Plot the Venn diagram
# plot(venn_data2)
plot_fig1b <- plot(venn_data2)

#=======================================================================

# ----- 5. Why summits matter (Figure 1c) -----

# Up to now:
# - You treated peaks as intervals

# Now:
#  - You treat peaks as points (summits)

# Why?
#  - Peak width varies wildly
#  - Binding precision lives at the summit
#  - Distance between summits is biologically meaningful

##### Figure 1c

# This figure requires a little more work. Let’s decompose it.

# Description of the figure in the paper:

# Position of TEAD4 peak summits relative to the summits of the overlapping YAP/TAZ peaks, 
# in a 500 bp window surrounding the summit of YAP/TAZ peaks.

# TAZ peaks coordinates and summit positions were used to represent common peaks 
# between YAP and TAZ peaks (YAP/TAZ peaks) and were used when comparing 
# YAP/TAZ peaks with other ChIP-seq data.

# - data we need: the TAZ peak set and the TEAD4 peak set.
# - x-axis: when a TEAD peak overlaps with a TAZ peak, 
#   the distance between the summit of the the TAZ peak and the TEAD4 peak summit.
# - y-axis: the number of TEAD peaks for each distance

# -------------------------------------------------------
# Each summit:
#  - Is a 1 bp interval
#  - Represents max ChIP signal

# A summit is the highest signal point within a peak. MACS3 outputs that.
#TAZ_summit<- import(here("data/fastq/TAZ_peak/TAZ_summits.bed"))
TAZ_summit<- import("N:/Bioinformatik/Programming/Projects/1_reproduce_figures_using_Supplementary_Data/data/fastq/TAZ_peak/TAZ_summits.bed")

# ----- 6. Restricting to biologically relevant summits -----
# Only keep TAZ summits from peaks that overlap YAP.
# So:
#  - You define YAP/TAZ peaks
#  - Using TAZ coordinates as the reference frame

# we will only subset the summit that is overlapping with YAP
## To make that line plot:
# 1. make sure two peaks first overlap with each other
# 2. calculate the summit distances between those two overlapping peaks
# TAZ_summit$name is <character> name
# name is the same name to subset the TAZ_summit
TAZ_summit<- TAZ_summit[
  TAZ_summit$name %in% TAZ_overlap_YAP_peaks$name
  ]

# TEAD4_summit<- import(here("data/fastq/TEAD4_peak/TEAD4_summits.bed"))
TEAD4_summit<- import("N:/Bioinformatik/Programming/Projects/1_reproduce_figures_using_Supplementary_Data/data/fastq/TEAD4_peak/TEAD4_summits.bed")

TEAD4_summit

# They both represent a single base point that has the highest signal in the peaks.
# -------------------------------------------------------

# ----- 7. Creating a ±250 bp window around TAZ summits -----
# This creates:
#  - A symmetric window
#  - From −250 bp to +250 bp
#  - Centered on TAZ summit

# expand the TAZ summit to a 500bp window
# because x-axis -250 to 250
# 250bp upstream and 250bp upstream
TAZ_500bp_window<- resize(TAZ_summit, width = 500, fix="center")

# ----- 8. Matching TEAD4 summits to TAZ windows -----
# Which TEAD4 summits fall into which TAZ windows
# One TEAD4 summit can match multiple TAZ windows
# One TAZ window can contain multiple TEAD4 summits
# findOverlaps of TEAD4_summit with TAZ_500bp_window
hits<- findOverlaps(TEAD4_summit, TAZ_500bp_window)

# a hits object with the indices of the overlapping query and subject
hits

# queryHits(hits) & subjectHits(hits) to extract the numbers
# calculating the distance between those two summits
# head(summit_distance)
summit_distance<- distance(TEAD4_summit[queryHits(hits)], TAZ_summit[subjectHits(hits)])

table(summit_distance)

# look at the distance which is equal to 0 & look at the TEAD4 summit
TEAD4_summit[queryHits(hits)][summit_distance ==0]

TAZ_summit[subjectHits(hits)][summit_distance ==0]

# The built-in `distance` function returns the pair-wise distance in absolute value.

# Let’s revise it to return negative values when TEAD4 summit precede the 
# TAZ summit and positive values when TEAD4 summit follows TAZ summit.


# ----- 9. Distance calculation (unsigned → signed) -----

## Unsigned distance
# summit_distance<- distance(TEAD4_summit[queryHits(hits)], TAZ_summit[subjectHits(hits)])

## Signed distance

# Interpretation:
#  - Negative → TEAD4 upstream of TAZ
#  - Positive → TEAD4 downstream of TAZ
#  - Zero → exact summit coincidence

# Compute signed distances
# create a new function
signed_distance <- function(A, B) {
  # Compute unsigned distance
  dist <- distance(A, B)
  
  # but if A is smaller than B
  # A is smaller than B and B-
  # If A is greater than B
  # Determine signs based on whether A precedes or follows B
  sign <- ifelse(start(A) < start(B), -1, 1)
  
  # Apply sign to distance
  dist * sign
}

library(dplyr)
library(ggplot2)
# summit_distance is a numeric vector
# Each element is a distance in base pairs
# Each number = one TEAD4–TAZ summit pair.
summit_distance<- signed_distance(TEAD4_summit[queryHits(hits)],
                                  TAZ_summit[subjectHits(hits)])
# table(summit_distance)
# Interpretation:
#  - Distance −249 bp occurred 1 time
#  - Distance -246 bp occurred 2 times
#  - Distance -233 bp occurred 1 time
#  - Distance -232 bp occurred 2 times

# So now you have:
#  - x-axis: distance
#  - y-axis: count (density proxy)

# table() returns a special object, not a normal data frame.

# class(table(summit_distance))

# A table:
#  - Looks like a matrix
#  - Behaves like a named vector
#  - Is awkward for ggplot2

# ggplot2 wants a data frame:
#  - one column per variable
#  - one row per observation

# So we need to convert!!!

# ----- 10. Turning distances into a density plot -----
# This converts:
#  - “many individual distances”
#     → “count per distance value”

# Then:
#  - Sort x-axis
# - Plot counts vs distance

# “Count how many summit pairs occur at each distance,
# then store the result in a tidy table so I can plot it.”

# This produces the raw spike plot, which is noisy.
# table() answers one question:
# “How many times does each value occur?”

# What %>% (the pipe) means
# “Take the result from the left and pass it as the first argument to the function on the right.”
distance_df<- table(summit_distance) %>%
  tibble::as_tibble() 
# So this is equivalent to:
  # tibble::as_tibble(table(summit_distance))

# What as_tibble() does
# tibble is a modern version of data.frame.
# tibble::as_tibble(table(summit_distance))
# Converts to:
# A tibble:
#   summit_distance     n
#   <chr>           <int>
# 1 -12                 2
# 2 0                   3
# 3 5                   3
# 4 10                  1

# Two columns:
#  - summit_distance → the distance (as character for now)
#  - n → how many times it occurred

# This format is exactly what ggplot() wants.

# Why summit_distance becomes character
# <summit_distance> <chr>
# Because:
#  - table() stores its names as strings
#  - as_tibble() preserves that

# That’s why later you had to do:
# mutate(summit_distance = as.numeric(summit_distance))

distance_df

#class(distance_df)
#str(distance_df)

# X-Axis | Y-Axis

# Plot it:
# It is exactly equivalent to:
# ggplot(distance_df, ...)

# It creates a plot object that says:
# “I intend to plot something using this data.
# The x-variable will be summit_distance.
# The y-variable will be n.”

distance_df %>%
  # aes() = aesthetic mapping
  # Map the column summit_distance to the x-axis
  # Map the column n to the y-axis
  ggplot(aes(x=summit_distance, y = n)) +
  geom_point()

# Hmm, something is off… the summit distance on the x-axis needs to be reordered

distance_df %>%
  # summit_distance (character because of table())
  # mutate() = add or transform columns in a data frame (or tibble).
  # Here, you’re overwriting summit_distance to be numeric.
  # class(distance_df$summit_distance)  # before: "character"
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  # Sorts the tibble by summit_distance ascending.
  # Result: rows go from the most negative distance to the most positive distance.
  arrange(summit_distance) %>%
  # class(distance_df$summit_distance)  # after mutate: "numeric"
  # Why important:
  #  - ggplot now treats x-axis as continuous instead of categorical.
  #  - The points will be positioned correctly along the genome, not evenly spaced like labels.
  
  # n (integer count of how many TEAD4–TAZ summit pairs are at each distance)
  # ggplot() creates a plot object.
  # aes() = aesthetic mapping (what goes where)
  # x = summit_distance → genomic distance relative to TAZ peak summit
  # y = n → number of TEAD4 summits at that distance
  ggplot(aes(x=summit_distance, y = n)) +
  # Actually draws a point for every row in the tibble:
  # x-coordinate = summit_distance (numeric now)
  # y-coordinate = n
  # Effectively: each dot = “TEAD4 summits count at this distance from TAZ peak.”
  geom_point()

# Let’s connect the points with a line:
# Take distance_df
distance_df %>%
  # Convert the distance column to numeric (so it’s a real axis)
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  # Sort the table so distances go from negative to positive
  arrange(summit_distance) %>%
  # Start a plot: x = distance, y = count
  ggplot(aes(x=summit_distance, y = n)) +
  # insdtead of geom_point() we use geom_line()
  # Draw one point per distance-count pair
  # This gives you a scatter plot of TEAD4 density around TAZ summits.
  geom_line()

# Result: super noisy / “wiggly” because counts jump a lot between neighboring base pairs

# ----- 11. Binning to smooth signal (5 bp bins) -----

# Groups distances into 5 bp bins
# Averages counts per bin
# Reduces noise
# Reveals structure

# This is effectively a discrete kernel smoothing without saying the word.

# Biologically:
#  - We don’t care about single-bp jitter
#  - We care about positional bias at ~10 bp scale

# The final plot answers:
# “Is TEAD4 centered on YAP/TAZ summits?”

# And the answer is: yes, very strongly.

# The plot looks too wigglely. Let’s smooth it by average the number of peaks per 5 bp bin.
df_binned <- distance_df %>%
  # Step 1: Convert distances to numeric & arrange
  mutate(summit_distance = as.numeric(summit_distance)) %>%
  arrange(summit_distance) %>%
  # summit_distance / 5
  # Divide each distance by 5 → groups of 5 bp (e.g., −12/5 = −2.4)
  # floor()
  # Round down to nearest integer (−2.4 → −3)
  # Ensures all distances in the same 5 bp interval map to the same bin
  # * 5
  # Convert back to original scale
  # −3 → −15, meaning bin covers −15 to −11
  # Result: a new column called bin, grouping distances into 5 bp windows:
  #   summit_distance	|   bin
  #         -12	      |   -15
  #         -11       | 	-15
  #         -10	      |   -10
  #         -9	      |   -10
  mutate(bin = floor(summit_distance / 5) * 5) %>%  # Create bins by grouping every 5 bp
  # Group by bin & summarize
  # group_by(bin) → treat all rows in the same 5 bp bin as one group
  group_by(bin) %>%
  # summarise(n = mean(n, na.rm=TRUE)) → compute average count in each bin
  # na.rm = TRUE ignores missing values
  # Result: one row per 5 bp bin, with smoothed y-values
  # summarise(n = mean(n, na.rm = TRUE))  # Calculate average 'n' for each bin
  summarise(n = mean(n, na.rm = TRUE)) %>%
  # Example:
  #    bin	|    n (mean)
  #   -250	|    0.6
  #   -245	|    1.2
  #   -240	|    1.0
  #    ...	|    ...
  # Use [[ ]] to force the column
  mutate(n = (n / max(n)) *2)   # scale 0–1
  # Effect: the wiggles are smoothed, trend is clear, small fluctuations are averaged out.

# Groups summit_distance values into 5 bp bins by dividing by 5, taking the floor (rounding down), and multiplying back by 5 to get the bin lower bound.
# View the binned dataframe
print(df_binned)

# 4: Plot the smoothed data
plot_fig1c <- df_binned %>%
  # x-axis = 5 bp bins
  # y-axis = average number of peaks in that bin
  ggplot(aes(x = bin, y = n)) +
  # Connect points with a line → shows trend of TEAD4 density around TAZ summits
  geom_line(color = "blue") +
  # Only label x-axis at −250, 0, +250 bp for readability
  scale_x_continuous(
    breaks = c(-250, 0, 250),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 2.125),
    breaks = seq(0, 2.125, 0.5),
    # sprintf("%.1f", x) → formats every number with 1 decimal place
    # So 0 → 0.0, 0.5 → 0.5
    labels = function(x) {
      sprintf("%.1f", x) %>%
        # sub("\\.0$", "", .) → removes the .0 at the end
        # So 0.0 becomes 0, 1.0 becomes 1
        sub("\\.0$", "", .)
    }, # remove .0 from whole numbers
    expand = c(0, 0) # start exactly at 0, no padding
  ) +
  # Human-readable axis labels
  xlab("Distance to the summit \nof TAZ peaks (bp)") +
  ylab("Peak density (x10²)") +
  # Clean background, larger font
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(fill = NA)) # keeps top/right borders


# Save each plot individually
# Figure 1a
png("fig_1_a.png", width = 10, height = 10, units = "in", res = 300)
plot(venn_data)
dev.off()

# Figure 1b
png("fig_1_b.png", width = 10, height = 10, units = "in", res = 300)
plot(venn_data2)
dev.off()

ggsave("fig_1_c.png", plot = plot_fig1c, width = 10, height = 10, units = "in")