```
Look at figure 1d, it is a scatter plot. what do we need?

we need the YAP, TAZ, and TEAD4 signal. This is the number of sequencing reads
x-axis: TEAD4 signal
y-axis: YAP1 signal
From the figure description:

> (d) Linear correlation between the signal of YAP or TAZ and TEAD4 peaks in the 5522 shared binding sites. r2 is the coefficients of determination of the two correlations.

Read in the peak files:
```{r}

library(GenomicRanges) # to deal with intervals
library(rtracklayer) # to read peak files or export peak files
library(here) # to construct the route path to read the files
library(dplyr)
library(ggplot2)

TAZ_peaks <- import(here("data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks <- import(here("data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peak<- import(here("data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

# find the overlaps
# it gives warnings, because there are unconventional chromosomes
YAP_overlap_TAZ_peaks <- subsetByOverlaps(YAP_peaks, TAZ_peaks)

# you can further overlap with TEAD4 peaks
YAP_overlap_TAZ_peaks_overlap_TEAD4 <- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)
YAP_overlap_TAZ_peaks_overlap_TEAD4

# you can export YAP_overlap_TAZ_peaks_overlap_TEAD4 GenomicRanges files into a .bed file
# use rtracklayer to write the GenomicRanges object to file
export(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
       con = here("data/fastq/YAP_TAZ_TEAD4_common.bed"))

# in wsl
# verify how many peaks are there
# wc -l YAP_TAZ_TEAD4_common.bed
## 5965 YAP_TAZ_TEAD4_common.bed

# look at the file
# head YAP_TAZ_TEAD4_common.bed
# 6 column bed file format
# .bed = min 3 columns: chromosome, start, end
## chr1    1437543 1438043 YAP_peak_1 100

# now you can actually count how many reads for each of this peak
# from those three transfer factor, actually .bam files

```
The next step is to get the ‘signal’ in those common peaks for YAP, TAZ and TEAD4, respectively. How do we do it?

The signal is the number of reads fall/mapped into those peaks/regions. and normalized to total number of reads (library size) for each experiment.

There are multiple ways to do it.

#### count the number of reads from bam files with bedtools
The [mutlicov](https://bedtools.readthedocs.io/en/latest/content/tools/multicov.html) subcommand from bedtools is what we need.

>bedtools multicov, reports the count of alignments from multiple position-sorted and indexed BAM files that overlap intervals in a BED file. Specifically, for each BED interval provided, it reports a separate count of overlapping alignments from each BAM file.

```{bash eval=FALSE}
cd data/fastq
# -bams = specify the bam files
# -bed = specify the intervals, common intervals, that we just overlapped
# then it will give you the Counts
bedtools multicov -bams YAP.sorted.bam TAZ.sorted.bam TEAD4.sorted.bam -bed YAP_TAZ_TEAD4_common.bed > YAP_TAZ_TEAD4_counts.tsv
```
It takes less than a minute to finish. Let’s take a look at the file
```{bash eval=FALSE}
head YAP_TAZ_TEAD4_counts.tsv
# chromosome  start   end     name    score   strand  YAP1_count  TAZ_count   TEAD4_count
# chr1    1024627 1025059 YAP_peak_3      494     .       88      72      82
# chr1    1264836 1265155 YAP_peak_4      148     .       32      37      88
# chr1    1265319 1265695 YAP_peak_5      131     .       29      31      26
# chr1    1360617 1360955 YAP_peak_6      306     .       46      52      88
# chr1    1659297 1659586 YAP_peak_8      45      .       15      14      20
# chr1    2061242 2061682 YAP_peak_10     356     .       54      65      60
# chr1    2140001 2140346 YAP_peak_11     86      .       27      18      27
# chr1    3543323 3543624 YAP_peak_12     155     .       24      30      28
# chr1    6724590 6724868 YAP_peak_14     251     .       38      42      90
# chr1    8061325 8061624 YAP_peak_17     62      .       21      34      38
```

# each experiment, or each transfer, they have different sequencing depths, 
# meaning they have different number of total reads in each experiment. 
# You want to normalize it to, for example, counts per million, 
# so you can kind of compare between different experiments. 

# we need to get the total number of reads in each  transferring factor chipseq

The last three columns are counts for YAP1, TAZ and TEAD4 in the common regions.

We need to normalize it to total number of reads in each library. Let’s use `samtools flagstat`:
```{bash eval=FALSE}

# Get all the numbers for those three bam files using same commands

# you can use samtools flagstat
# conda env list = environment list
# $ samtools flagstat YAP.sorted.bam # in wsl
# Output:
# total number of reads: 24549590
# 24549590 + 0 in total (QC-passed reads + QC-failed reads)
# 24549590 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# mapped reads:
# we need to use this number to normalize the counts
# take YAP_count and divide it by 24549590, then multiply by 1 million to get CPM
# so we get counts per million (CPM) normalized counts
# 23653961 + 0 mapped (96.35% : N/A)
# 23653961 + 0 primary mapped (96.35% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)


# $ samtools flagstat TAZ.sorted.bam
# 27521260 + 0 in total (QC-passed reads + QC-failed reads)
# 27521260 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 26789648 + 0 mapped (97.34% : N/A)
# 26789648 + 0 primary mapped (97.34% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

# $ samtools flagstat TEAD4.sorted.bam
# 34776462 + 0 in total (QC-passed reads + QC-failed reads)
# 34776462 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 34332907 + 0 mapped (98.72% : N/A)
# 34332907 + 0 primary mapped (98.72% : N/A)
# 0 + 0 paired in sequencing
# 0 + 0 read1
# 0 + 0 read2
# 0 + 0 properly paired (N/A : N/A)
# 0 + 0 with itself and mate mapped
# 0 + 0 singletons (N/A : N/A)
# 0 + 0 with mate mapped to a different chr
# 0 + 0 with mate mapped to a different chr (mapQ>=5)

```
So the total number of priamry mapped reads are: 23653961, 26789648 and 34332907 for YAP, TAZ and TEAD4, respectively.

Load the data into R:
```{r}

# install.packages("readr")
library(readr)
# total number of mapped reads
counts<- read_tsv(here("data/fastq/YAP_TAZ_TEAD4_counts.tsv"), col_names = FALSE)
colnames(counts)<- c("chr", "start", "end", "name", "score", "value", "YAP1", "TAZ", "TEAD4")

head(counts)
# > head(counts)
# # A tibble: 6 × 9
#   chr     start     end name        score value  YAP1   TAZ TEAD4
#   <chr>   <dbl>   <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl>
# 1 chr1  1024627 1025059 YAP_peak_3    494 .        88    72    82
# 2 chr1  1264836 1265155 YAP_peak_4    148 .        32    37    88
# 3 chr1  1265319 1265695 YAP_peak_5    131 .        29    31    26
# 4 chr1  1360617 1360955 YAP_peak_6    306 .        46    52    88
# 5 chr1  1659297 1659586 YAP_peak_8     45 .        15    14    20
# 6 chr1  2061242 2061682 YAP_peak_10   356 .        54    65    60

```
normalize the counts to `CPM` (counts per million)
```{r}

counts<- counts %>%
  mutate(YAP1 = YAP1/23653961 * 10^6,
         TAZ = TAZ/26789648 * 10^6,
         TEAD4 = TEAD4/34332907 * 10^6)

head(counts)
```
ready to plot
```{r}

# # A tibble: 6 × 9
#   chr     start     end name        score value  YAP1   TAZ TEAD4
#   <chr>   <dbl>   <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl>
# 1 chr1  1024627 1025059 YAP_peak_3    494 .     3.72  2.69  2.39
# 2 chr1  1264836 1265155 YAP_peak_4    148 .     1.35  1.38  2.56
# 3 chr1  1265319 1265695 YAP_peak_5    131 .     1.23  1.16  0.757
# 4 chr1  1360617 1360955 YAP_peak_6    306 .     1.94  1.94  2.56
# 5 chr1  1659297 1659586 YAP_peak_8     45 .     0.634 0.523 0.583
# 6 chr1  2061242 2061682 YAP_peak_10   356 .     2.28  2.43  1.75

# YAP1 vs TEAD4 column
ggplot(counts, aes(x = TEAD4, y = YAP1)) +
    geom_point()
# there is one outlier, one big signal in top right corner

counts %>%
    filter(TEAD4 > 60)
#        this is the peak reagion
#              v       v
#   chr      start      end name          score value  YAP1   TAZ TEAD4
#   <chr>    <dbl>    <dbl> <chr>         <dbl> <chr> <dbl> <dbl> <dbl>
# 1 chr14 61754882 61755803 YAP_peak_3030  4298 .      82.6  80.3  64.5

```
Whenever you see something abnormal, always be suspisious.
Try to look at IGV, wether there is something strange there.
There is an outlier with strong signal (note, check it on IGV to see if it is real, it could be a black-listed region with strong signal)

in IGV it looks real and is not one of the blacklist regions from ENCODE.

NOTE: Download the blacklisted regions from here: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz

We can remove that outlier, or use log2 scale

```{r}

# log2 scale
ggplot(counts, aes(x = TEAD4, y = YAP1)) +
    geom_point(color = "#ff4000") +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    theme_classic(base_size = 40) +
    # add some labels
    xlab("TEAD4 signal") +
    ylab("YAP1 signal")

# ggplot(counts, aes(x = TEAD4, y = YAP1)) +
#     geom_point(color = "#ff4000") +
#     scale_x_continuous(trans = "log2") +
#     scale_y_continuous(trans = "log2") +
#     theme_classic(base_size = 14) +
#     # add some labels
#     xlab("TEAD4 signal") +
#     ylab("YAP1 signal")

```
We still need is a regretion line and r squared value

We will use [ggpmisc](https://github.com/aphalo/ggpmisc) to add the R^2.

```{r}

# install.packages("ggpmisc")
# library(ggpmisc)
# ggplot(counts, aes(x=TEAD4, y= YAP1)) +
#   geom_point(color = "#ff4000") +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
#   # this function is from ggpmisc package to add the R^2
#   stat_poly_eq(
#     aes(label = ..rr.label..),
#     formula = y ~ x,
#     parse = TRUE,
#     color = "black"
#   ) +
#   scale_x_continuous(trans = 'log2') +
#   scale_y_continuous(trans = 'log2') +
#   theme_classic(base_size = 14) +
#   xlab("TEAD4 signal") +
#   ylab("YAP1 signal")

library(ggpmisc)
ggplot(counts, aes(x = TEAD4, y = YAP1)) +
    geom_point(color = "#ff4000") +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") + # Linear regression line
    # this function is from ggpmisc package to add the R^2
    stat_poly_eq(
        aes(label = ..rr.label..),
        formula = y ~ x,
        parse = TRUE,
        color = "black",
        size = 12
    ) +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 35) +
    xlab("TEAD4 signal") +
    ylab("YAP signal")

# ggsave("fig_1_d_1.png", width = 10, height = 10, units = "in")

```

correlation coefficent is the r which ranges from -1 to 1. Coefficient of Determination is the R^2.

```{r}

# calculate correlation
correlation_coefficent<- cor(log2(counts$TEAD4), log2(counts$YAP1))
correlation_coefficent
# 0.8095894

# square it
R_squared<- correlation_coefficent^2
R_squared
# 0.655435
```
We can re-create the other scatter plot:
```{r}
ggplot(counts, aes(x = TEAD4, y = TAZ)) +
    geom_point(color = "#ff4000") +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") + # Linear regression line
    stat_poly_eq(
        aes(label = ..rr.label..),
        formula = y ~ x,
        parse = TRUE,
        color = "black",
        size = 12
    ) +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    coord_fixed(ratio = 1) +
    theme_classic(base_size = 35) +
    xlab("TEAD4 signal") +
    ylab("TAZ signal")

# ggsave("fig_1_d_2.png", width = 10, height = 10, units = "in")

# ggplot(counts, aes(x=TEAD4, y= TAZ)) +
#   geom_point(color = "#ff4000") +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +  # Linear regression line
#   stat_poly_eq(
#     aes(label = ..rr.label..),
#     formula = y ~ x,
#     parse = TRUE,
#     color = "black"
#   ) +
#   scale_x_continuous(trans = 'log2') +
#   scale_y_continuous(trans = 'log2') +
#   theme_classic(base_size = 14) +
#   xlab("TEAD4 signal") +
#   ylab("TAZ signal")
```
Tip: take a look at [ggpubr](https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html)
```{r}



### Figure 1 f

```
>f. Absolute distance of YAP peaks (n=7709), TAZ peaks (n=9798), TEAD4 peaks (n=8406) or overlapping YAP/TAZ/TEAD peaks (n=5522) to the nearest TSS.

Figure 1f is a stacked bar plot. It shows the proportion of the peaks grouped by their distance to the closest TSS (transcription start site).

We need to get the distance of YAP/TAZ/TEAD4 peaks to the nearest TSS.
And then get the table.
Make that Stacked Barplot.

I will show you how to do this from scratch:

```{r}

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(GenomicFeatures)

# Get all the transcripts from database:
# TxDb.Hsapiens.UCSC.hg38.knownGene
# Get the TSS
hg38_transcripts <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38_transcripts

# transcript =  11121-14413
# transcript name = ENST00000832824.1

# We also need to get the transcription start site (TSS)
# get the TSS.
tss_gr <- promoters(hg38_transcripts, upstream = 0, downstream = 1)

# first one is 11121
# last one is 324923
# GRanges object with 412044 ranges and 2 metadata columns:
#                       seqnames    ranges strand |     tx_id           tx_name
#                          <Rle> <IRanges>  <Rle> | <integer>       <character>
#        [1]                chr1     11121      + |         1 ENST00000832824.1
#        [2]                chr1     11125      + |         2 ENST00000832825.1
#        [3]                chr1     11410      + |         3 ENST00000832826.1
#        [4]                chr1     11411      + |         4 ENST00000832827.1
#        [5]                chr1     11426      + |         5 ENST00000832828.1
#        ...                 ...       ...    ... .       ...               ...
#   [412040] chrX_MU273397v1_alt    260095      - |    412040 ENST00000710260.1
#   [412041] chrX_MU273397v1_alt    282686      - |    412041 ENST00000710028.1
#   [412042] chrX_MU273397v1_alt    316302      - |    412042 ENST00000710030.1
#   [412043] chrX_MU273397v1_alt    315236      - |    412043 ENST00000710216.1
#   [412044] chrX_MU273397v1_alt    324923      - |    412044 ENST00000710031.1

# Calculate the distance to the nearest TSS
# distanceToNearest(GenomicRanges object, GenomicRanges object with len=1 because 1 base TSS)
distance_to_tss <- distanceToNearest(YAP_peaks, tss_gr)

# Print the distance
distance_to_tss

# Hits object 
# Hits object with 9807 hits and 1 metadata column:
#  "index YAP_peaks" "index tss_gr"    "distance"
#          queryHits subjectHits |  distance
#          <integer>   <integer> | <integer>
#      [1]         1           1 |       796
#      [2]         2          46 |       125
#      [3]         3         349 |      4504
#      [4]         4       18448 |      1267
#      [5]         5       18448 |      1750
#      ...       ...         ... .       ...
#   [9803]      9803      378030 |      4939
#   [9804]      9804      383719 |     18165
#   [9805]      9805      383732 |      2477
#   [9806]      9806      378136 |      6869
#   [9807]      9807      383802 |      1165
```
It is a Hits object, and we can access the the distance metadata column
```{r}

# access distance column
mcols(distance_to_tss)
# DataFrame with 9807 rows and 1 column
#       distance
#      <integer>
# 1          796
# 2          125
# 3         4504
# 4         1267
# 5         1750
# ...        ...
# 9803      4939
# 9804     18165
# 9805      2477
# 9806      6869
# 9807      1165

# you get a vector of distances
head(mcols(distance_to_tss)$distance)
# [1]  796  125 4504 1267 1750  821

```
Let’s do that for all three factors:
```{r}

YAP_dist<- mcols(distanceToNearest(YAP_peaks, tss_gr))$distance
TAZ_dist<- mcols(distanceToNearest(TAZ_peaks, tss_gr))$distance
TEAD4_dist<- mcols(distanceToNearest(TEAD4_peak, tss_gr))$distance

```
put them in a single dataframe
```{r}

# put it in the same dataframe
# convert to data.frame first (data.frame(...))
# data.frame(factor = "YAP", distance = YAP_dist)
tss_distance_df<- bind_rows(data.frame(factor = "YAP", distance = YAP_dist),
          data.frame(factor = "TAZ", distance = TAZ_dist),
          data.frame(factor = "TEAD4", distance = TEAD4_dist))
          
head(tss_distance_df)
#   factor distance
# 1    YAP      796
# 2    YAP      125
# 3    YAP     4504
# 4    YAP     1267
# 5    YAP     1750
# 6    YAP      821

# create the categories
tss_distance_df %>%
    mutate(category = case_when(
        distance < 1000 ~ "<1kb",
        distance >= 1000 & distance < 10000 ~ "1-10kb",
        distance >= 10000 & distance <= 100000 ~ "10-100kb",
        distance > 100000 ~ "100kb"
    )) %>%
    head()

#   factor distance category
# 1    YAP      796     <1kb
# 2    YAP      125     <1kb
# 3    YAP     4504   1-10kb
# 4    YAP     1267   1-10kb
# 5    YAP     1750   1-10kb
# 6    YAP      821     <1kb

```
You can see how I build the pipe %>% step by step.
```{r}

# group by factor and category, then count
counts_per_category<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  group_by(factor, category) %>%
  count()

counts_per_category

# A tibble: 12 × 3
# Groups:   factor, category [12]
#    factor category     n
#    <chr>  <chr>    <int>
#  1 TAZ    1-10kb    4423
#  2 TAZ    10-100kb  3632
#  3 TAZ    <1kb      2556
#  4 TAZ    >100kb     108
#  5 TEAD4  1-10kb    4856
#  6 TEAD4  10-100kb  4029
#  7 TEAD4  <1kb      2490
#  8 TEAD4  >100kb     137
#  9 YAP    1-10kb    3929
# 10 YAP    10-100kb  3197
# 11 YAP    <1kb      2596
# 12 YAP    >100kb      85


# Original figure is percentage
# You need the total number of peaks for each category
total_counts<- tss_distance_df %>%
  mutate(category = case_when(
    distance < 1000 ~ "<1kb",
    distance >=1000 & distance < 10000 ~ "1-10kb",
    distance >= 10000 & distance <=100000 ~ "10-100kb",
    distance > 100000 ~ ">100kb"
  )) %>%
  count(factor, name = "total")

total_counts

#   factor total
# 1    TAZ 10719
# 2  TEAD4 11512
# 3    YAP  9807

```
You need to merge those two
```{r}
# join = merge
merged_df<- left_join(counts_per_category, total_counts)
merged_df %>%
    # calculate percentage
    mutate(Percentage = n / total * 100) %>%
    # plot using ggplot2
    ggplot(aes(x = factor, y = Percentage, fill = category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Distance to TSS",
        x = "Group",
        y = "Percentage"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_classic(base_size = 14)

```
You can customize the color and reorder the category as you want.
```{r}
merged_df$category<- factor(merged_df$category, 
                            levels = c("<1kb", "1-10kb", "10-100kb", ">100kb"))
merged_df %>%
    mutate(Percentage = n / total * 100) %>%
    ggplot(aes(x = factor, y = Percentage, fill = category)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        title = "Distance to TSS",
        x = "Group",
        y = "Percentage"
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(values = c("#EF3E2B", "#F16161", "#F59595", "#FAD1C8")) +
    theme_classic(base_size = 14)

# ggsave("fig_1_f.png", width = 10, height = 10, units = "in")