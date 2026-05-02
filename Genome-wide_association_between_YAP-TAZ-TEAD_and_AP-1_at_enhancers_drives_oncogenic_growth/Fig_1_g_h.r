```
1. Analyse what kind of plot type are they.

    g) IGV track view for raw signal file. BigWig file (just a histogram).
    h) Pie Chart for all TAZ/YAP/TEAD4 peaks.
       Annotate them, whether they are:
       - promoters
       - active enhancers
       - inactive enhancers
       - not classified
```{r}

```
Figure 1g

Figure 1g is an IGV track view. We already generated the bigwig files and we can easily load them to IGV and take a screenshot.

Below is one of the known YAP1 target gene CCN2: IGV_view.png

The paper also shows H3K4me1, H3K4me3 and H3K27ac tracks which the author analyzed from a [previous dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49651) for MDA-MB-231 breast cancer cell line.

I will leave the analysis from the raw fastq files as an exercise for you.
```{r}

```
plot the tracks from scratch

If you look carefully, the genome track view for the bigwig files is just a histogram.

If you want to plot the bigwig signal tracks:

take a look at [plotgardner](https://phanstiellab.github.io/plotgardener/articles/guides/plotting_multiomic_data.html) bioconductor package.

[karyoploteR](https://bioconductor.org/packages/release/bioc/html/karyoploteR.html) I used for my [scATACutils package](https://github.com/crazyhottommy/scATACutils) for plotting [scATACseq tracks](https://rpubs.com/crazyhottommy/scATAC_tracks).

[Gviz](https://bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html)

[ggbio](https://www.bioconductor.org/packages/release/bioc/html/ggbio.html)

[pyGenometracks](https://github.com/deeptools/pyGenomeTracks) in python.
```{r}

```
Figure 1h

Figure 1h is a pie chart showing the proportion of the annotation of the YAP/TAZ/TEAD4 common peaks. The annotations are: promoters, active enhancers, inactive enhancers and not classified.

promoters are YAP/TAZ/TEAD4 peaks overlapping with H3K4me3 peaks
active enhancers are YAP/TAZ/TEAD4 peaks overlapping with H3K4me1 and H3K27ac peaks
inactive enhancers are YAP/TAZ/TEAD4 peaks overlapping with H3K4me1 but not H3K37ac peaks
Read this classic ChIP-seq paper from Keji Zhao’s group to understand what genomic features those histone modifications are associated with.

To do this, I will download the H3K4me1, H3K4me3 and H3K27ac peaks from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49651.

click the ftp or http in the link above, and you will download the peaks.txt.gz files.

Note, the data were aligned to hg19 human genome. We did our alignment using hg38 genome. We will need to lift-over the hg19 peak files to hg38 version.

ls data/public_data

GSE49651_MDAMB231.H3K27Ac.hg19.all.MACS.1e-5_peaks.txt.gz GSE49651_MDAMB231.H3K4me3.hg19.ALL.MACS.1e-5_peaks.txt.gz
GSE49651_MDAMB231.H3K4me1.hg19.all.MACS.1e-5_peaks.txt.gz

Those are the .bed files of the peaks called by MACS2.

Look at one of the files:
$ zless -S GSE49651_MDAMB231.H3K27Ac.hg19.all.MACS.1e-5_peaks.txt.gz
chr     start   end     name            score
chr1    713349  715152  MACS_peak_1     280.72
chr1    753156  753893  MACS_peak_2     80.69
chr1    761553  763199  MACS_peak_3     110.65
chr1    839067  841109  MACS_peak_4     399.13
chr1    893458  897246  MACS_peak_5     885.66
chr1    900528  903402  MACS_peak_6     3048.62
chr1    910760  912393  MACS_peak_7     56.71
chr1    933792  938040  MACS_peak_8     1735.46
chr1    947481  951255  MACS_peak_9     2086.04
chr1    954216  961837  MACS_peak_10    3100.00
chr1    993622  995996  MACS_peak_11    396.06
chr1    996781  1000389 MACS_peak_12    491.03
chr1    1003863 1005560 MACS_peak_13    96.73
chr1    1014575 1016235 MACS_peak_14    1454.98
chr1    1050277 1052815 MACS_peak_15    747.73
chr1    1092421 1094966 MACS_peak_16    3100.00
chr1    1165532 1168829 MACS_peak_17    1090.38


PROBLEM:

aligned this old data sets using hg38 human genome, the peak files are in hg19 version.
We need to lift-over the hg19 peak files to hg38 version.


To do this, we can use the liftOver tool from UCSC genome browser.

Download the chain file from UCSC 
https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/ and 
the command line tool from https://hgdownload.soe.ucsc.edu/admin/exe/

You need to also download the chain file which maps those two genomes: hg19 and hg38.
So you know the region hg19 maps to which region in hg38.
```{r}

```
curl -O https://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# I am using a macs with arm64 system
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/liftOver
curl -O https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
# make it executable
chmod u+x liftOver
```{r}

```
Now, we are ready to lift-over the coordinates from hg19 to hg38:

liftOver   input.bed                                                 chainfile                output.bed  unmapped.bed
./liftOver GSE49651_MDAMB231.H3K27Ac.hg19.all.MACS.1e-5_peaks.txt.gz hg19ToHg38.over.chain.gz H3K27ac.bed unmapped1.bed

./liftOver GSE49651_MDAMB231.H3K4me3.hg19.ALL.MACS.1e-5_peaks.txt.gz hg19ToHg38.over.chain.gz H3K4me3.bed unmapped2.bed

./liftOver GSE49651_MDAMB231.H3K4me1.hg19.all.MACS.1e-5_peaks.txt.gz hg19ToHg38.over.chain.gz H3K4me1.bed unmapped3.bed

# how many peaks after the liftover, lines in each bed file
 wc -l *bed
   37480 H3K27ac.bed
   26765 H3K4me1.bed
   17991 H3K4me3.bed
      30 unmapped1.bed
      56 unmapped2.bed
      12 unmapped3.bed
```{r}

```
Again, we can use bedtools to do all those intersections and annotate the YAP1/TAZ/TEAD4 peaks. I will keep everything in R instead:

```{r}

library(rtracklayer) # will import those files as GenomeRanges objects
library(here)
library(dplyr)
library(ggplot2)

TAZ_peaks<- import(here("data/fastq/TAZ_peak/TAZ_peaks.narrowPeak"))
YAP_peaks<- import(here("data/fastq/YAP_peak/YAP_peaks.narrowPeak"))
TEAD4_peak<- import(here("data/fastq/TEAD4_peak/TEAD4_peaks.narrowPeak"))

# overlap with histone modification peaks
YAP_overlap_TAZ_peaks<- subsetByOverlaps(YAP_peaks, TAZ_peaks)

YAP_overlap_TAZ_peaks_overlap_TEAD4<- subsetByOverlaps(YAP_overlap_TAZ_peaks, TEAD4_peak)
YAP_overlap_TAZ_peaks_overlap_TEAD4

#> GRanges object with 5965 ranges and 6 metadata columns:
#>          seqnames              ranges strand |          name     score
#>             <Rle>           <IRanges>  <Rle> |   <character> <numeric>
#>      [1]     chr1     1024628-1025059      * |    YAP_peak_3       494
#>      [2]     chr1     1264837-1265155      * |    YAP_peak_4       148
#>      [3]     chr1     1265320-1265695      * |    YAP_peak_5       131
#>      [4]     chr1     1360618-1360955      * |    YAP_peak_6       306
#>      [5]     chr1     1659298-1659586      * |    YAP_peak_8        45
#>      ...      ...                 ...    ... .           ...       ...
#>   [5961]     chrX 154368850-154369243      * | YAP_peak_9801        90
#>   [5962]     chrX 154596614-154596846      * | YAP_peak_9802       120
#>   [5963]     chrX 154600351-154600918      * | YAP_peak_9803       131
#>   [5964]     chrX 154732680-154732891      * | YAP_peak_9804        58
#>   [5965]     chrX 155888248-155888487      * | YAP_peak_9806       108
#>          signalValue    pValue    qValue      peak
#>            <numeric> <numeric> <numeric> <integer>
#>      [1]    16.37320   54.3883  49.46070       192
#>      [2]     8.81841   18.8472  14.89750       188
#>      [3]     8.62413   17.0083  13.14010       112
#>      [4]    14.37360   35.1297  30.64260       148
#>      [5]     5.06808    7.8940   4.58709       102
#>      ...         ...       ...       ...       ...
#>   [5961]     6.13997  12.64820   9.00333       237
#>   [5962]     8.21346  15.84150  12.02720       128
#>   [5963]     8.62413  17.00830  13.14010       417
#>   [5964]     5.74942   9.32421   5.89513        78
#>   [5965]     6.72473  14.58900  10.84320       151
#>   -------
#>   seqinfo: 27 sequences from an unspecified genome; no seqlengths

```
read in the histone modification peaks:
```{r}

H3K4me1<- import(here("data/public_data/H3K4me1.bed"))
H3K4me3<- import(here("data/public_data/H3K4me3.bed"))
H3K27ac<- import(here("data/public_data/H3K27ac.bed"))

```
From the method section: How the authors defined where there is an enhancer, active enhancer or promoter:

The presence of H3K4me1 and H3K4me3 peaks, their genomic locations and their overlap were the 
criteria used to define promoters and enhancers: 
i) H3K4me3 peaks not overlapping with H3K4me1 peaks and close to a TSS (± 5kb) were defined as 
promoters, as NA otherwise; ii) H3K4me1 peaks not overlapping with H3K4me3 peaks were defined 
as enhancers; iii) regions with the co-presence of H3K4me1 and H3K4me3 peaks were visually 
inspected on IGV and were defined as promoters, enhancers or NA after the evaluation of the 
proximity to a TSS and the comparison of the enrichment signals. Finally, promoters or enhancers 
were defined as active if overlapping with H3K27ac peaks.

Define (in)active enhancers:
```{r}

# subsetByOverlaps() is a GenomicRanges function (Bioconductor) that 
# filters one set of genomic ranges based on whether they overlap 
# another set of genomic ranges.
# “Keep (or remove) ranges from x depending on whether they overlap ranges in y.”
# subsetByOverlaps(x, y, invert = FALSE)
# x: GRanges object to be filtered
# y: GRanges object used as the reference
# invert = FALSE (default):
# keep ranges in x that overlap y
# invert = TRUE:
# keep ranges in x that do NOT overlap y
active_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac)
inactive_enhancers<- subsetByOverlaps(H3K4me1, H3K27ac, invert=TRUE)

promoters<- subsetByOverlaps(H3K4me3, H3K4me1, invert=TRUE)

```
YAP/TAZ/TEAD peaks annotation YAP/TAZ/TEAD peaks were annotated as promoters or enhancers if 
their summit was overlapping with promoter or enhancer regions as defined above. 
Peaks with the summit falling in regions with no H3K4me1 or H3K4me3 peaks, or in NA 
regions were defined as “not assigned” and discarded from subsequent analyses.
```{r}

n_active_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                      active_enhancers) %>%
  length()

n_inactive_enhancers<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4,
                                        inactive_enhancers) %>%
  length()

n_promoters<- subsetByOverlaps(YAP_overlap_TAZ_peaks_overlap_TEAD4, 
                               promoters) %>%
  length()

n_unclassified<- length(YAP_overlap_TAZ_peaks_overlap_TEAD4) - n_active_enhancers -
  n_inactive_enhancers - n_promoters

```
Put the numbers in a dataframe:
```{r}

annotation_df<- data.frame(category = c("active_enhancers", "inactive_enhancers",
                        "promoters", "unclassified"),
peak_number = c(n_active_enhancers, n_inactive_enhancers, 
                n_promoters, n_unclassified))


annotation_df

```
Make the pie chart:
```{r}
library(ggplot2)

ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  # category section is a Barplot
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_brewer(palette = "Set3") 

```
change the order of the categories by changing the factor level
```{r}

annotation_df$category<- factor(annotation_df$category, 
                                levels = c("promoters", "active_enhancers",
                                           "inactive_enhancers", "unclassified"))

colors<- c("#8D1E0F", "#F57D2B", "#FADAC4", "#D4DADA")

ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors)

```
Note, the authors re-called the peaks using their own IgG sample:

Peak calls and read density tracks were generated using SPP version 1.1148 
with default parameters and using as control sample the IgG ChIP-seq data 
generated in our laboratory because of the low sequencing depth of 
the Input DNA contained in SRP028597.

This can cause drastically different number of H3K4me1/3 and H3K27ac peaks. 
It is not surprising to me that we now have a lot of unclassified peaks.
```{r}

# Add the percentage to the pie chart
# Calculate percentages and cumulative positions for labeling
annotation_df <- annotation_df %>%
  dplyr::mutate(
    percentage = peak_number / sum(peak_number) * 100,
    label = paste0(round(percentage, 1), "%")
  )

annotation_df

# Create the pie chart
ggplot(annotation_df, aes(x = "", y = peak_number, fill = category)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar("y", start = 0) +
  theme_void() + # Remove unnecessary axes
  labs(title = "YAP/TAZ/TEAD4 peaks") +
  scale_fill_manual(values = colors) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 7) # Add percentage labels

ggsave("fig_1_h.png", width = 10, height = 10, units = "in", bg = "white")