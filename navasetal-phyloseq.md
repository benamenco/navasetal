
<link href="http://joey711.github.com/phyloseq/markdown.css" rel="stylesheet"></link>

# Support for QIIME in R


## Why use R

First published in 1996, ["R"](http://www.R-project.org) is an integrated software application and programming language designed for interactive data analysis (R Core Team). It is available for Linux, Mac OS, and Windows free of charge under an open-source license [GPL2](http://www.gnu.org/licenses/gpl-2.0.html). Since its inception R has had an emphasis on being a tool for interactive statistical analysis through [functional programming](http://en.wikipedia.org/wiki/Functional_programming), where the primary investigative and inferrential activites are performed by writing a series of repeatable commands into "scripts" that can be recorded and published. This paradigm lends itself well to reproducible research, and is enhanced substantially by R's integration with tools for literate programming such as [Sweave](http://www.stat.uni-muenchen.de/~leisch/Sweave/) ( Leisch 2002), [knitr](http://cran.r-project.org/web/packages/knitr/citation.html) (Xie 2013), and R [markdown](http://cran.r-project.org/web/packages/markdown/index.html) (Allaire 2013), as well as data graphics. There are thousands of free and open-source extensions to R ("packages") available from the main R repository, [CRAN](http://cran.r-project.org/), further organized by volunteer experts into 31 [task views](http://cran.r-project.org/web/views/). Among these are dedicated package views relevant to microbiome data, including [phylogenetics](http://cran.r-project.org/web/views/Phylogenetics.html), [clustering](http://cran.r-project.org/web/views/Cluster.html), [environmetrics](http://cran.r-project.org/web/views/Environmetrics.html), [machine learning](http://cran.r-project.org/web/views/MachineLearning.html), [multivariate](http://cran.r-project.org/web/views/Multivariate.html) and [spatial statistics](http://cran.r-project.org/web/views/Spatial.html), as well as a separate reviewed and curated repository dedicated to biological statistics called [Bioconductor](http://bioconductor.org/) (over 600 packages).
	
	
## Support for QIIME output files in phyloseq 

At present, support for QIIME in R is predominantly achieved through a package called ["phyloseq"](http://dx.plos.org/10.1371/journal.pone.0061217) (McMurdie and Holmes 2013), dedicated to the reproducible analysis of microbiome census data in R. phyloseq defines an object-oriented data class for the consistent representation of related (heterogenous) microbiome census data that is independent of the sequencing or OTU-clustering method (storing OTU abundance, taxonomy classification, phylogenetic relationships, representative biological sequences and sample covariates). The package supports QIIME by including functions for importing data from biom-format files derived from more recent versions of QIIME (`import_biom`) as well as legacy OTU-taxonomy delimited files (`import_qiime` and related user accessible subfunctions). Later editions of phyloseq (>`1.5.15`) also include an API for importing data directly from the [microbio.me/qiime](http://www.microbio.me/qiime/) data repository. In all cases, these API functions return an instance of the "phyloseq" class that contains the available heterogenous components in "native" R classes. phyloseq includes a number of tools for connecting with other microbiome analysis functions available in other R packages, as well as its own functions for flexible graphics production built using [ggplot2](http://ggplot2.org/) (Wickham 2009), demonstrated in supplemental files (File SXX) and [online tutorials](http://joey711.github.io/phyloseq/). For researchers interested in developing or using methods not directly supported by phyloseq, nor its data infrastructure, the biom-format specific core functions in phyloseq have been migrated to an official API in the biom-format project as an installable R package called ["biom"](http://cran.r-project.org/web/packages/biom/), now released on CRAN. This also includes some biom-format specific functionality that is beyond the scope of phyloseq, though support for QIIME is still likely best achieved using phyloseq. 


---

##  Examples

As with some of the earlier examples of QIIME commands with corresponding output and figures, in this section we have included some key R commands potentially useful during interactive analysis in the R environment. For simplicity we are only showing results related to the open-reference OTU data, stored in an object in our examples named `open`, and imported into R using the phyloseq command `import_biom`.


```r
open = import_biom("path-to-file.biom", ...)
```


Additional input data files can also be provided to `import_biom`, or merged with `open` after its instantiation. For clarity, subsets and transformations of the data in `open` are stored in objects having names that begin with `open`. As with the remainder of the examples highlighted in this section, the complete code sufficient for reproducing all results and figures are included in the R Markdown originated document, Figure SNN, which also includes several additional examples and details not shown here.

Although not always very illuminating, a comparison of OTU-richness between samples or groups of samples can easily be achieved with the `plot_richness` command. For the most precise estimates of richness for most samples, this should be performed *before* random subsampling or other transformations of the abundance data. Here `open` contains data that has already been randomly subsampled. In the following graphic we can see that the wild type samples are generally more diverse (higher richness) and somewhat more variable than the transgenic samples for essentially all body sites, though the differences between the two mice genotypes are small.


```r
plot_richness(open, x = "BODY_SITE", color = "GENOTYPE") + geom_boxplot()
```


This plot command also illustrates the use of a function in ggplot2, `geom_boxplot`, that instructs the ggplot2 graphics engine to add an additional graphical element -- in this case a boxplot for each of the natural groups in the graphic. These available additional graphical instructions (called "layers" in the grammar of graphics nomenclature) are embedded with the returned plot object for subsequent rendering, inspection, or further modification, allowing for powerfully customized representations of the data.

Here is an example leveraging the abundance bar plot function from phyloseq, `plot_bar`, in order to compare the relative abundances of key phyla between the wild type and transgenic mice across body sites. The first step was actually some additional data transformations (not shown, see Figure SNN) in order to subset the data to only major expected phyla (`subset_taxa`), merge OTUs from the same phyla as one entry (`merge_taxa`), and merge samples from the same body site and mouse genotype (`merge_samples`).

In this first version fo the figure, a base graphic is created in which the bars are color-filled according to phyla, with each rectangle representing the relative abundance of a phylum in a particular sample group and then stacked from largest to smallest. The `facet_grid(~GENOTYPE)` layer then separates the data into two equivalent panels, one for each genotype in the data.  


```r
p2 = plot_bar(openphyab, "bodysite", fill = "phyla", title = title)
p2 + facet_grid(~GENOTYPE)
```


From this first bar plot it is clear that all body sites from the average wildtype mouse have Firmicutes as their phylum of largest cumulative proportion, except for the "feces", where it is anyway a close call between firmicutes and bacteroidetes. By contrast, some of the average transgenic mice samples have a much higher proportion of proteobacteria or bacteroidetes than the corresponding wild type samples. One draw back of this kind of stacked bar representation is that it is difficult to compare any of the sub-bars except for those at the bottom. This can be easily alleviated by changing the `facet_grid` call such that a separate panel is made for each phyla in the dataset, as follows. 


```r
p2 + facet_grid(phyla ~ GENOTYPE) + ylim(0, 100)
```


With essentially the same effort to produce, the 14 panels of this second bar plot graphic allow an easy and quantitative comparison of the relative abundances of each phyla across body sites and genotype.


#### Ordination methods.

Microbiome datasets can be quite multivariate in nature, and dimensional reduction ([ordination](http://en.wikipedia.org/wiki/Ordination)) methods can be a useful form of exploratory analysis to better understand some of the largest patterns in the data. Many ordination methods are wrapped in phyloseq by the `ordinate` function, and many more are offered in available R packages. Here is an example performing multidimensional scaling (MDS) on the precomputed unweighted UniFrac distance matrix for the open-reference dataset. The ordination result (`openUUFMDS`) is first passed to `plot_scree` in order to explore the "scree plot" representing the relative proportions of variability represented by each successive axis. Both the ordination result and the original data are then passed to `plot_ordination` with sufficent parameters to shade the sample points by genotype, and create separate panels for each body site. 


```r
openUUFMDS = ordinate(open, "MDS", distance = UniFrac[["unweighted"]][["open"]])
plot_scree(openUUFMDS, "Eignevalues for unweighted UniFrac MDS")
plot_ordination(open, openUUFMDS, color = "GENOTYPE") + geom_point(size = 5) + 
    facet_wrap(~BODY_SITE)
```





It appears that a subset of the wildtype samples from all but the mouth and abdomen-skin body sites cluster toward the left of the plot. This appears to be the major pattern along the axis that also comprises the greatest proportion of variability in the dataset. At this stage of analysis it seems worthwhile to try to identify which OTU abundances are most different between these groups, and then perform some formal validation/testing of these differences.


#### Heatmaps

Heatmaps can be useful for exploratory analysis of microbiomes by mapping abundance values to a color scale in a condensed, pattern-rich format. A good heatmap graphic will often generate hyotheses about sample and/or OTU clustering in the data that must then be more fully vetted. There are two structural aspects of a heatmap graphic that control whether or not it will reveal any interpretable patterns: (1) the ordering of axes, and (2) the color scaling. For example, the `plot_heatmap()` command in the phyloseq package takes a similar approach to NeatMap (Rajaram and Oono, PMID 20096121), in that it uses ordination results rather than hierarchical clustering to determine the index order for each axis. For `plot_heatmap`, the default color scaling maps a particular shade of blue to a log transformation of abundance that works well for microbiome data, although this transformation is also user-definable. 

In creating this example heatmap, a key step was proper filtering of the data. We removed OTUs that do not appear in very many samples. The possible contribution to the graphic of these seldom-occurring OTUs is limited, more often contributing to "noise" that causes the heatmap to look dark, empty, and uninterpretable (see example in File SNN). We used a non-metric multidimensional scaling of the Bray-Curtis distance to determine the order of the OTUs and samples. From this representation, it is possible to distinguish high-level patterns and simultaneously note the samples and OTUs involved. For instance, all but a few of the mouth samples are in a cluster toward the right of the heatmap. One of the key features of this group is an obvious relative overabundance of three firmicutes OTUs that are among some of the most abundant in this subset of the data. Similarly, another clear pattern is a distinction between a group of wild type samples form various body sites on the left of the heatmap that appear to have higher proportions of a number of different firmicutes OTUs, as well as a few specific bacteroidetes OTUs. This is distinct from -- and an approximate reversal of -- the largest cluster of samples in the center-right of the graphic, in which many of the most-abundant OTUs are a different subset of bacteroidetes OTUs and just a few firmicutes. We also found it helpful to further pursue these high-level patterns by splitting the data into firmicutes-only and bacteroidetes-only subsets, and then plotting new heatmaps with finer-scale taxonomic labels. This required essentially the same commands and limited additional effort, well tailored for exploratory interactive analysis, much of which we have documented in File SNN.

The call to create Figure NNN was the following in R/phyloseq, omitting some details to improve the axis labels for publication:


```r
title = "plot_heatmap using NMDS/Bray-Curtis for both axes ordering"
plot_heatmap(openfpp, "NMDS", "bray", taxa.label = "Phylum", sample.label = "bsgt", 
    title = title)
```



Caption
**Figure NNN. Example heatmap of the high-level patterns in the open-reference dataset**. The graphic was produced by the `plot_heatmap()` function in phyloseq after subsetting the data to the most-prevalent 100 OTUs (See File SNN). The order of sample and OTU elements was determined by the radial position of samples/OTUs in the first two axes of a Non-metric multidimensional scaling (NMDS) of the Bray-Curtis distance. Other choices for distance and ordination method can also be useful. The horizontal axis represents samples, with the genotype and body site labeled, while the vertical axis represents OTUs, labeled by phyla. Both axes are further color-coded to emphasize the different categories of labels. The blue-shade color scale indicates the abundance of each OTU in each sample, from black (zero, not observed) to very light blue (highly abundant, >1000 reads).



---

## Begin Evaluated Code Examples

This document serves as a reproducible, supplementary tutorial for importing the example data and performing some example analyses in R using phyloseq.

The following includes some examples importing QIIME data and creating data graphics. Some preprocessing steps are omitted because they were performed prior to import using tools in QIIME.

Use the `import_biom` function in phyloseq to import `.biom` and associated files from each dataset type into R, represented as a single integrated object. First we have to load the phyloseq extension/library/package (these terms are roughly equivalent in R).

Load required packages for import, show package version.

```r
library("phyloseq")
packageVersion("phyloseq")
```

```
## [1] '1.5.20'
```



### Define files

Find the data files that we want to import. Uses some regular expression matching recursively in the directory we're in.


```r
biomfiles = list.files(getwd(), "\\.biom$", recursive = TRUE)
names(biomfiles) <- c("closed", "denovo", "open")
# Find the phylogenetic tree files.
treefiles = list.files(getwd(), "\\.tre[e]{0,}$", recursive = TRUE)
names(treefiles) <- c("denovo", "closed", "open")
# Find the DNA sequence (`.fasta`) file(s) with the representative
# sequences for each OTU.  Only available for 'closed'
rsfiles = list.files(getwd(), "\\.fasta$", recursive = TRUE)
names(rsfiles) <- c("closed")
# Find and import the 'mapping file' with sample covariate data.  There's
# only one.
mapfiles = list.files(getwd(), "\\_map\\.txt$", recursive = TRUE)
sada = import_qiime_sample_data(mapfiles)
```



### Import data

In this case, we show as an example importing the open-reference files.


```r
openrefRfile = "imported-open-ref.RData"
```


If the files have not already been imported and saved yet then do that. If the result file is already present in compressed `.RData` format, then this import part is skipped. To re-run, you can simply delete the file, imported-open-ref.RData.


```r
if (!file.exists(openrefRfile)) {
    open = import_biom(biomfiles["open"], treefiles["open"], parseFunction = parse_taxonomy_greengenes)
    # Add sample map data to the open-reference data object called `open`.
    open = merge_phyloseq(open, sada)
    # Remove any samples that don't have a body site `NA`.
    open = subset_samples(open, !is.na(BODY_SITE))
    # Save to openrefRfile
    save(open, file = openrefRfile)
} else {
    load(openrefRfile)
}
```



---

### Import closed-reference or de novo files

One could repeat the process above to import the files for closed-reference and de novo OTU clustering, already identified in the file name vectors above. However, these comparisons are already described in other parts of the article to which this is a supplement, and do not need to be repeated here.


---

### Import pre-calculated UniFrac distance matrices


```r
UUFfiles = list.files(getwd(), "^unweighted\\_unifrac\\_dm\\.txt$", recursive = TRUE)
WUFfiles = list.files(getwd(), "^weighted\\_unifrac\\_dm\\.txt$", recursive = TRUE)
names(UUFfiles) <- c("closed", "denovo", "open")
names(WUFfiles) <- c("closed", "denovo", "open")
# Read each distance matrix as a data.frame, store in a list
UUFL = lapply(UUFfiles, read.table, header = TRUE)
WUFL = lapply(WUFfiles, read.table, header = TRUE)
# Re-order the open-reference distances to match data object `open`
reorderdists = function(distdf, newI) {
    distdf <- distdf[newI, ]
    distdf <- distdf[, newI]
    return(distdf)
}
UUFL = lapply(UUFL, reorderdists, newI = sample_names(open))
WUFL = lapply(WUFL, reorderdists, newI = sample_names(open))
# Coerce each data.frame to a 'dist' object
UUFL = lapply(UUFL, as.dist)
WUFL = lapply(WUFL, as.dist)
# Combine into a bigger nested list
UniFrac = list(unweighted = UUFL, weighted = WUFL)
```



---

### load ggplot2 and do themes

Load the ggplot2 package for full access of its graphics tools and customizations. Then set some custom default plot parameters for graphics theme and discrete color scale.


```r
library("ggplot2")
theme_set(theme_bw())
scale_colour_discrete <- function(...) {
    scale_colour_brewer(palette = "Set1", ...)
}
scale_fill_discrete <- function(...) {
    scale_fill_brewer(palette = "Set1", ...)
}
```



### Richness estimates

Create a summary of richness estimates, attempting different combinations of variable organization. Simply plotting all of the samples, even with labels, would not be very illuminating. Here we try a few different choices for organization, and find that it looks like the microbiomes associated with the WT genotype have a richness that is more variable and possibly more diverse (larger).


```r
# Remove the 'UBERON:' label from each BODY_SITE entry
sample_data(open)$BODY_SITE <- gsub("^UBERON\\:", "", get_variable(open, "BODY_SITE"))
plot_richness(open, x = "BODY_PRODUCT") + geom_boxplot()
```

![plot of chunk richness-comparisons](figure/richness-comparisons1.png) 

```r
prich = plot_richness(open, x = "BODY_SITE", color = "GENOTYPE")
prich
```

![plot of chunk richness-comparisons](figure/richness-comparisons2.png) 

```r
# Additional ggplot2 adjustments for publication
prich = prich + geom_boxplot()
prich = prich + theme(text = element_text(size = 12), axis.text.x = element_text(vjust = 0.5))
prich = prich + ggtitle("OTU Richness of Samples, Grouped by Body Site")
prich
```

![plot of chunk richness-comparisons](figure/richness-comparisons3.png) 

```r
# Command to save the figure as PDF, then JPG
dev.off()
```

```
## null device 
##           1
```

```r
ggsave("Figure_22.pdf", prich, width = 7.5, height = 5.25)
ggsave("Figure_22.jpg", prich, width = 7.5, height = 5.25, dpi = 300)
```


### Flexible plots of abundances
I want to make a plot of relative abundance by `BODY_SITE`, but I also want to be able to show the effect of `GENOTYPE`, being informed by some earlier exploratory analysis which is show later. The first thing to do is to add a dummy variable for merging that combines the two categorical variables I want to compare.


```r
# Make a new lower-case bodysite label with 'UBERON:' removed
sample_data(open)$bodysite <- gsub("^UBERON\\:", "", get_variable(open, "BODY_SITE"))
# change 'multi-tissue structure' to just multi-tissue
sample_data(open)$bodysite <- factor(gsub(" structure$", "", get_variable(open, 
    "bodysite")))
# Merge based on the dummy variable, bsgt
sample_data(open)$bsgt <- paste0(get_variable(open, "bodysite"), get_variable(open, 
    "GENOTYPE"))
sample_data(open)$bsgt <- factor(sample_data(open)$bsgt)
openbsgtmerge = merge_samples(open, "bsgt")
```

```
## Warning: NAs introduced by coercion
```


Repair R factors that were damaged (coerced to integers) during standard merge.


```r
# repair factors
sample_data(openbsgtmerge)$bodysite <- levels(sample_data(open)$bodysite)[get_variable(openbsgtmerge, 
    "bodysite")]
sample_data(openbsgtmerge)$GENOTYPE <- levels(sample_data(open)$GENOTYPE)[get_variable(openbsgtmerge, 
    "GENOTYPE")]
# transform merged data to proportions out of 100 (percent)
openbsgtmerge = transform_sample_counts(openbsgtmerge, function(x) 100 * x/sum(x))
```


Define a special set of Phyla that we want to focus on in our plotting. We have prior knowledge that these are important, and that we might want to lump the other phyla into an "other" category.


```r
plotphyla = c("Proteobacteria", "Fusobacteria", "Firmicutes", "Deferribacteres", 
    "Cyanobacteria", "Bacteroidetes", "Actinobacteria")
tax_table(openbsgtmerge) <- cbind(phyla = NA_character_, tax_table(openbsgtmerge))
wh1 = which(tax_table(openbsgtmerge)[, "Phylum"] %in% plotphyla)
tax_table(openbsgtmerge)[wh1, "phyla"] <- tax_table(openbsgtmerge)[wh1, "Phylum"]
```


Now use `plot_bar` to show the relative abundances. We can show with and without separate by `GENOTYPE`.

```r
title = "OTU abundance; phylum, body site, and genotype"
p1 = plot_bar(openbsgtmerge, "bodysite", fill = "phyla", title = title)
p1 + facet_grid(~GENOTYPE)
```

![plot of chunk plotbar0](figure/plotbar01.png) 

```r
p1 + facet_grid(phyla ~ GENOTYPE) + ylim(0, 100)
```

![plot of chunk plotbar0](figure/plotbar02.png) 


Now the same style of plot, but after removing the unknown (NA) OTUs, and agglomerating the OTUs from the same phylum. Figure 23.


```r
title = "Phylum abundance; body site, and genotype"
openphyab = openbsgtmerge
tax_table(openphyab) <- tax_table(openphyab)[, c("Kingdom", "phyla")]
openphyab = tax_glom(openphyab, "phyla")
p2 = plot_bar(openphyab, "bodysite", fill = "phyla", title = title)
p2 = p2 + theme(text = element_text(size = 12), axis.text.x = element_text(vjust = 0.5))
# Figure 23
fig23 = p2 + facet_grid(~GENOTYPE)
fig23
```

![plot of chunk Figure-23](figure/Figure-23.png) 

```r
dev.off()
```

```
## null device 
##           1
```

```r
ggsave("Figure_23.pdf", fig23, width = 7.5, height = 5.25)
ggsave("Figure_23.jpg", fig23, width = 7.5, height = 5.25, dpi = 300)
```


Figure 24 - faceting. Created by adding layers to previous figure, `p2`.


```r
fig24 = p2 + facet_grid(phyla ~ GENOTYPE, scales = "free_y")
fig24 = fig24 + theme(axis.text.y = element_text(size = 7), strip.text.y = element_text(size = 6))
fig24
```

![plot of chunk Figure-24](figure/Figure-24.png) 

```r
dev.off()
```

```
## null device 
##           1
```

```r
ggsave("Figure_24.pdf", fig24, width = 7.5, height = 7)
ggsave("Figure_24.jpg", fig24, width = 7.5, height = 7, dpi = 300)
```


Can also try to visualize on a phylogenetic tree, since we have one integrated with this data (including all the transformations, filtering we've been doing)


```r
plot_tree(openphyab, nodelabf = nodeplotboot(), color = "bodysite", label.tips = "phyla", 
    shape = "GENOTYPE", size = "abundance", ladderize = "left", plot.margin = 0.45)
```

![plot of chunk plot-trees](figure/plot-trees.png) 



### MDS/PCoA on weighted and unweighted UniFrac distances

In phyloseq, many native dimensional reduction ([ordination](http://en.wikipedia.org/wiki/Ordination)) methods are wrapped by the `ordinate` function, which calculates and returns different ordination results classes, depending on the `method` chosen for ordination (see `help("ordinate")` for more details).

Here we can visually compare the MDS ordinations on weighted and unweighted UniFrac distances.


```r
openUUFMDS = ordinate(open, "MDS", distance = UniFrac[["unweighted"]][["open"]])
openWUFMDS = ordinate(open, "MDS", distance = UniFrac[["weighted"]][["open"]])
plot_scree(openUUFMDS, "Unweighted UniFrac MDS")
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open1.png) 

```r
plot_scree(openWUFMDS, "weighted UniFrac MDS")
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open2.png) 

```r
plot_ordination(open, openUUFMDS, color = "BODY_SITE") + geom_point(size = 5)
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open3.png) 

```r
pu = plot_ordination(open, openUUFMDS, color = "GENOTYPE") + geom_point(size = 5)
pu
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open4.png) 

```r
pu + facet_wrap(~BODY_SITE)
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open5.png) 

```r
plot_ordination(open, openWUFMDS, color = "BODY_SITE") + geom_point(size = 5)
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open6.png) 

```r
pw = plot_ordination(open, openWUFMDS, color = "GENOTYPE") + geom_point(size = 5)
pw
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open7.png) 

```r
pw + facet_wrap(~BODY_SITE)
```

![plot of chunk plot-ordination-open](figure/plot-ordination-open8.png) 


Grab the left-hand subset of data, label cluster.


```r
pus = pu
pus$data <- pu$data[which(pu$data$Axis.1 < -0.05 & pu$data$Axis.2 > -0.25), 
    ]
pus + facet_wrap(~BODY_SITE)
```

![plot of chunk subset-UUFMDS-testing](figure/subset-UUFMDS-testing.png) 

```r
# Add 'cluster' label to sample data
sample_data(open)$cluster <- "right"
sample_data(open)[rownames(pus$data), "cluster"] <- "left"
```


Testing taxa abundance correlation between binary groups using permutation t-test and multiple hypothesis correction, `mt`. For speed and power, we first apply a filter in which OTUs are removed that do not change much across samples, standard deviation less than 5 (arbitrary, not rigorously chosen threshold).


```r
# hist(apply(otu_table(open), 1, sd), breaks=500, xlim=c(10, 200),
# ylim=c(0, 50))
testOTUs = names(which(apply(otu_table(open), 1, sd) > 5))
openfilt = prune_taxa(testOTUs, open)
opmt = mt(openfilt, "cluster")
```

```
## B=10000
## b=100	b=200	b=300	b=400	b=500	b=600	b=700	b=800	b=900	b=1000	
## b=1100	b=1200	b=1300	b=1400	b=1500	b=1600	b=1700	b=1800	b=1900	b=2000	
## b=2100	b=2200	b=2300	b=2400	b=2500	b=2600	b=2700	b=2800	b=2900	b=3000	
## b=3100	b=3200	b=3300	b=3400	b=3500	b=3600	b=3700	b=3800	b=3900	b=4000	
## b=4100	b=4200	b=4300	b=4400	b=4500	b=4600	b=4700	b=4800	b=4900	b=5000	
## b=5100	b=5200	b=5300	b=5400	b=5500	b=5600	b=5700	b=5800	b=5900	b=6000	
## b=6100	b=6200	b=6300	b=6400	b=6500	b=6600	b=6700	b=6800	b=6900	b=7000	
## b=7100	b=7200	b=7300	b=7400	b=7500	b=7600	b=7700	b=7800	b=7900	b=8000	
## b=8100	b=8200	b=8300	b=8400	b=8500	b=8600	b=8700	b=8800	b=8900	b=9000	
## b=9100	b=9200	b=9300	b=9400	b=9500	b=9600	b=9700	b=9800	b=9900	b=10000	
## r=1	r=2	r=3	r=4	r=5	r=6	r=7	r=8	r=9	r=10	
## r=11	r=12	r=13	r=14	r=15	r=16	r=17	r=18	r=19	r=20	
## r=21	r=22	r=23	r=24	r=25	r=26	r=27	r=28	r=29	r=30	
## r=31	r=32	r=33	r=34	r=35	r=36	r=37	r=38	r=39	r=40	
## r=41	r=42	r=43	r=44	r=45	r=46	r=47	r=48	r=49	r=50	
## r=51	r=52	r=53	r=54	r=55	r=56	r=57	r=58	r=59	r=60	
## r=61	r=62	r=63	r=64	r=65	r=66	r=67	r=68	r=69	r=70	
## r=71	r=72	r=73	r=74	r=75	r=76	r=77	r=78	r=79	r=80	
## r=81	r=82	r=83	r=84	r=85	r=86	r=87	r=88	r=89	r=90	
## r=91	r=92	r=93	r=94	r=95	r=96	r=97	r=98	r=99	r=100	
## r=101	r=102	r=103	r=104	r=105	
```

```r
opmt[opmt$adjp < 0.01, c("adjp", "Phylum", "Class", "Genus")]
```

```
##                       adjp         Phylum                 Class
## 308076              0.0028  Bacteroidetes           Bacteroidia
## 431900              0.0028  Bacteroidetes           Bacteroidia
## 212201              0.0028  Bacteroidetes           Bacteroidia
## 199118              0.0028  Bacteroidetes           Bacteroidia
## 272850              0.0028  Bacteroidetes           Bacteroidia
## 355746              0.0028  Bacteroidetes           Bacteroidia
## 261177              0.0028  Bacteroidetes           Bacteroidia
## New.ReferenceOTU182 0.0028           <NA>                  <NA>
## 374370              0.0028  Bacteroidetes           Bacteroidia
## 277624              0.0028  Bacteroidetes           Bacteroidia
## 398839              0.0028 Proteobacteria    Betaproteobacteria
## New.ReferenceOTU61  0.0028  Bacteroidetes           Bacteroidia
## 262938              0.0028     Firmicutes       Erysipelotrichi
## 181719              0.0028  Bacteroidetes           Bacteroidia
## 228329              0.0028  Bacteroidetes           Bacteroidia
## 187028              0.0028  Bacteroidetes           Bacteroidia
## 174056              0.0028  Bacteroidetes           Bacteroidia
## New.ReferenceOTU176 0.0028           <NA>                  <NA>
## New.ReferenceOTU189 0.0028     Firmicutes            Clostridia
## 325918              0.0028     Firmicutes               Bacilli
## 339593              0.0028     Firmicutes               Bacilli
## New.ReferenceOTU232 0.0028  Bacteroidetes           Bacteroidia
## 340960              0.0028     Firmicutes               Bacilli
## 260582              0.0028 Proteobacteria   Deltaproteobacteria
## New.ReferenceOTU195 0.0028  Bacteroidetes           Bacteroidia
## 266060              0.0028  Bacteroidetes           Bacteroidia
## 393399              0.0028  Bacteroidetes           Bacteroidia
## New.ReferenceOTU169 0.0028 Proteobacteria                  <NA>
## 192222              0.0028  Bacteroidetes           Bacteroidia
## 231755              0.0028  Bacteroidetes           Bacteroidia
## 1857626             0.0028     Firmicutes            Clostridia
## 261908              0.0028     Firmicutes            Clostridia
## New.ReferenceOTU138 0.0028     Firmicutes            Clostridia
## 586001              0.0028    Tenericutes            Mollicutes
## 2729098             0.0028 Proteobacteria Epsilonproteobacteria
## 180671              0.0028     Firmicutes            Clostridia
## 239562              0.0028     Firmicutes            Clostridia
## New.ReferenceOTU47  0.0028     Firmicutes            Clostridia
## 268733              0.0028     Firmicutes            Clostridia
## 177352              0.0028     Firmicutes            Clostridia
## 701462              0.0028     Firmicutes            Clostridia
## 3998912             0.0028     Firmicutes            Clostridia
## 383716              0.0028           <NA>                  <NA>
## 1107439             0.0060     Firmicutes            Clostridia
##                             Genus
## 308076                       <NA>
## 431900                       <NA>
## 212201                       <NA>
## 199118                       <NA>
## 272850                       <NA>
## 355746                       <NA>
## 261177                       <NA>
## New.ReferenceOTU182          <NA>
## 374370                       <NA>
## 277624                 Prevotella
## 398839                 Sutterella
## New.ReferenceOTU61           <NA>
## 262938                Allobaculum
## 181719                Bacteroides
## 228329                       <NA>
## 187028                       <NA>
## 174056                       <NA>
## New.ReferenceOTU176          <NA>
## New.ReferenceOTU189          <NA>
## 325918              Lactobacillus
## 339593              Lactobacillus
## New.ReferenceOTU232    Prevotella
## 340960              Lactobacillus
## 260582              Desulfovibrio
## New.ReferenceOTU195          <NA>
## 266060                       <NA>
## 393399                       <NA>
## New.ReferenceOTU169          <NA>
## 192222                 Prevotella
## 231755                       <NA>
## 1857626                      <NA>
## 261908                       <NA>
## New.ReferenceOTU138          <NA>
## 586001                       <NA>
## 2729098              Helicobacter
## 180671                       <NA>
## 239562                       <NA>
## New.ReferenceOTU47           <NA>
## 268733                Clostridium
## 177352                       <NA>
## 701462                      rc4-4
## 3998912                      <NA>
## 383716                       <NA>
## 1107439                      <NA>
```


Permanova / adonis


```r
library("vegan")
adonis(UniFrac$unweighted$open ~ cluster + GENOTYPE + BODY_SITE, as(sample_data(open), 
    "data.frame"))
```

```
## 
## Call:
## adonis(formula = UniFrac$unweighted$open ~ cluster + GENOTYPE +      BODY_SITE, data = as(sample_data(open), "data.frame")) 
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)    
## cluster    1      2.41   2.405   13.85 0.139  0.001 ***
## GENOTYPE   1      0.22   0.215    1.24 0.012  0.153    
## BODY_SITE  6      2.85   0.475    2.74 0.165  0.001 ***
## Residuals 68     11.81   0.174         0.683           
## Total     76     17.28                 1.000           
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


Mantel test

```r
mantel(UniFrac$unweighted$open, UniFrac$weighted$open)
```

```
## 
## Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel(xdis = UniFrac$unweighted$open, ydis = UniFrac$weighted$open) 
## 
## Mantel statistic r: 0.688 
##       Significance: 0.001 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0751 0.1032 0.1250 0.1496 
## 
## Based on 999 permutations
```



### Heatmap
The `plot_heatmap` function in phyloseq can create heatmaps that are very useful for exploratory analysis, leading you to hypotheses that you can then more fully vet. For more details and examples, see [the plot_heatmap tutorial](http://joey711.github.io/phyloseq/plot_heatmap-examples). Briefly, there are two extremely important aspects of a heatmap graphic that can often dictate whether or not it reveals any interpretable patterns: (1) the ordering of axes, and (2) the color scaling. Addressing item 1, the `plot_heatmap` function in phyloseq takes a similar approach to [the NeatMap package](http://cran.r-project.org/web/packages/NeatMap/index.html), in that it uses ordination results rather than hierarchical clustering to determine the axis orders. Addressing item 2, the color scaling maps color shading to a particular log transformation of abundance that seems to work well for microbiome data, and is also completely user-definable with all manner of alternative scalings supported by [the scales package](http://cran.r-project.org/web/packages/scales/index.html), including linear and identity transformations. However, we find that some manner of log transformation is usually best.

To make an interpretable heatmap, we need to do some filtering of OTUs that just don't appear in very many samples and make the graphic look dark and empty. They are not contributing much to our understanding of the samples, so why should we include them in the graphic? Here is what I mean (for fairness, using the `plot_heatmap` function and MDS/UniFrac as the ordination method):


```r
plot_heatmap(open, "MDS", UniFrac$unweighted$open, sample.label = "GENOTYPE")
```

![plot of chunk heatmap0](figure/heatmap0.png) 


By default `plot_heatmap` omits the labels of axes with a large number of features. Even on a screen with fairly high resolution it is difficult to discern any but a few of the most abundant and prevalent OTUs in the previous graphic.

Let's instead trim OTUs that do not appear in very many samples in the dataset. We'll sort first by prevalance (defined here as number of samples appeared) and then by abundance value. This new filtered version of the dataset is called `openfp`.

We need to decide how many OTUs to include, and this should be based on the data in some way. For example here, I've plotted the prevalence of OTUs to help justify/indicate why 100 OTUs was chosen. It is an arbitrary value on its own, but conservative when considering the number of samples in which OTUs appear (see plot below). In a sorted order by prevalence, the 100th OTU is already appearing in less than 26% of samples. Adding more OTUs is adding sparsity to a graphic in which sparsity makes patterns more difficult to see. The following is the code that filters these 100 "most prevalent" OTUs, and stores this filtered data as `openfp`.


```r
samobs = apply(otu_table(open), 1, function(x, m) sum(x > m), m = 2)
main = "Number of samples in which OTU was observed more than twice"
plot(sort(samobs, TRUE)[1:350], type = "h", main = main, ylab = "Number of samples", 
    xlab = "OTU order")
```

![plot of chunk trim-by-prevalence](figure/trim-by-prevalence.png) 

```r
otudf = data.frame(prev = samobs, sums = taxa_sums(open))
otudf = otudf[order(-otudf$prev, -otudf$sums), ]
openfp = prune_taxa(rownames(otudf)[1:100], open)
```


Count and show total reads for each phylum...

```r
tapply(taxa_sums(openfp), factor(tax_table(openfp)[, "Phylum"]), sum)
```

```
##  Actinobacteria   Bacteroidetes   Cyanobacteria Deferribacteres 
##            1799           30846             221             867 
##      Firmicutes  Proteobacteria     Tenericutes             TM7 
##           30766           10940             198             118
```

```r
get_taxa_unique(openfp, "Phylum")
```

```
## [1] "Bacteroidetes"   "Proteobacteria"  "Firmicutes"      "Cyanobacteria"  
## [5] NA                "Tenericutes"     "Actinobacteria"  "TM7"            
## [9] "Deferribacteres"
```


Looks like most of the remaining reads/OTUs in this dataset are Bacteroidetes, Firmicutes, Proteobacteria, and Actinobacteria. Subset to those. Also, change the combined label from `bsgt` to `body site and genotype`.

```r
keepPhy = c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria")
openfpp = subset_taxa(openfp, Phylum %in% keepPhy)
openfpp
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 93 taxa and 77 samples ]
## sample_data() Sample Data:       [ 77 samples by 56 sample variables ]
## tax_table()   Taxonomy Table:    [ 93 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 93 tips and 92 internal nodes ]
```


Yep. Only lost 7 OTUs.

#### plot_heatmap

Now let's use this filtered data to create a more interpretable heatmap, still using `plot_heatmap` in phyloseq.


```r
title = "plot_heatmap using MDS/UniFrac for sample axis order"
plot_heatmap(openfpp, "MDS", UniFrac$unweighted$open, taxa.label = "Phylum", 
    sample.label = "GENOTYPE", title = title)
```

![plot of chunk explore-heatmaps-p3](figure/explore-heatmaps-p31.png) 

```r
title = "plot_heatmap using NMDS/UniFrac for sample axis order"
plot_heatmap(openfpp, "NMDS", UniFrac$unweighted$open, taxa.label = "Phylum", 
    sample.label = "GENOTYPE", title = title)
```

![plot of chunk explore-heatmaps-p3](figure/explore-heatmaps-p32.png) 

```r
title = "plot_heatmap using NMDS/Bray-Curtis for both axes ordering"
p3 = plot_heatmap(openfpp, "NMDS", "bray", taxa.label = "Phylum", sample.label = "GENOTYPE", 
    title = title)
p3
```

![plot of chunk explore-heatmaps-p3](figure/explore-heatmaps-p33.png) 

```r
p3 + facet_grid(~GENOTYPE)
```

![plot of chunk explore-heatmaps-p3](figure/explore-heatmaps-p34.png) 

```r
p3 + facet_grid(Phylum ~ GENOTYPE)
```

![plot of chunk explore-heatmaps-p3](figure/explore-heatmaps-p35.png) 

```r
p3 + facet_grid(Phylum ~ .)
```

![plot of chunk explore-heatmaps-p3](figure/explore-heatmaps-p36.png) 


For publication, re-create the combined version, and also subset each phyla and re-make the corresponding individual heatmaps for each of Bacteroidetes and Proteobacteria. 


```r
xlab = "body site and genotype"
title = "plot_heatmap using NMDS/Bray-Curtis for both axes ordering"
p4 = plot_heatmap(openfpp, "NMDS", "bray", taxa.label = "Phylum", sample.label = "bsgt", 
    title = title)
p4$scales$scales[[1]]$name <- xlab
p4
```

![plot of chunk article-heatmaps-p4](figure/article-heatmaps-p41.png) 

```r
openfpp.bac = subset_taxa(openfpp, Phylum == "Bacteroidetes")
title = "plot_heatmap, NMDS/Bray-Curtis, Bacteroidetes only"
p4 = plot_heatmap(openfpp.bac, "NMDS", "bray", taxa.label = "Family", sample.label = "bsgt", 
    title = title)
p4$scales$scales[[1]]$name <- xlab
p4
```

![plot of chunk article-heatmaps-p4](figure/article-heatmaps-p42.png) 

```r
openfpp.fir = subset_taxa(openfpp, Phylum == "Firmicutes")
title = "plot_heatmap, NMDS/Bray-Curtis, Firmicutes only"
p4 = plot_heatmap(openfpp.fir, "NMDS", "bray", taxa.label = "Family", sample.label = "bsgt", 
    title = title)
p4$scales$scales[[1]]$name <- xlab
p4
```

![plot of chunk article-heatmaps-p4](figure/article-heatmaps-p43.png) 


Colorize axis labels? The ggplot2 package doesn't yet support an aesthetic mapping between a data variable and text color. Here is the added code for doing that (and results).


```r
# Full data
title = "plot_heatmap using NMDS/Bray-Curtis for both axes ordering"
p6 = plot_heatmap(openfpp, "NMDS", "bray", taxa.label = "Phylum", sample.label = "bsgt", 
    title = title)
p6 = p6 + xlab("samples: body site and genotype")
gt = gsub("^([[:print:]]+)(WT|TG)$", "\\1", p6$scales$scales[[1]]$labels)
xfac = factor(gt)
p7 = p6 + theme(axis.text.x = element_text(colour = (scales::brewer_pal("qual", 
    1))(nlevels(xfac))[xfac]))
yfac = factor(p6$scales$scales[[2]]$labels)
p7 = p7 + theme(axis.text.y = element_text(colour = (scales::brewer_pal("qual", 
    2))(nlevels(yfac))[yfac]))
p7$scales$scales[[1]]$name <- xlab
p7 + theme(text = element_text(size = 16), axis.text.x = element_text(size = 12), 
    axis.text.y = element_text(size = 12))
```

![plot of chunk colorize-axis-labels](figure/colorize-axis-labels1.png) 

```r
# Bacteroidetes only
title = "plot_heatmap, NMDS/Bray-Curtis, Bacteroidetes only"
p6 = plot_heatmap(openfpp.bac, "NMDS", "bray", taxa.label = "Family", sample.label = "bsgt", 
    title = title)
p6 = p6 + xlab("samples: body site and genotype")
gt = gsub("^([[:print:]]+)(WT|TG)$", "\\2", p6$scales$scales[[1]]$labels)
xfac = factor(gt)
p7 = p6 + theme(axis.text.x = element_text(colour = (scales::brewer_pal("qual", 
    2))(nlevels(xfac))[xfac]))
yfac = factor(p6$scales$scales[[2]]$labels)
p7 = p7 + theme(axis.text.y = element_text(colour = (scales::brewer_pal("qual", 
    3))(nlevels(yfac))[yfac]))
p7$scales$scales[[1]]$name <- xlab
p7
```

![plot of chunk colorize-axis-labels](figure/colorize-axis-labels2.png) 

```r
# Firmicutes only
title = "plot_heatmap, NMDS/Bray-Curtis, Firmicutes only"
p6 = plot_heatmap(openfpp.fir, "NMDS", "bray", taxa.label = "Family", sample.label = "bsgt", 
    title = title)
p6 = p6 + xlab("samples: body site and genotype")
gt = gsub("^([[:print:]]+)(WT|TG)$", "\\2", p6$scales$scales[[1]]$labels)
xfac = factor(gt)
p7 = p6 + theme(axis.text.x = element_text(colour = (scales::brewer_pal("qual", 
    2))(nlevels(xfac))[xfac]))
yfac = factor(p6$scales$scales[[2]]$labels)
p7 = p7 + theme(axis.text.y = element_text(colour = rep((scales::brewer_pal("qual", 
    1))(8), 2)[1:nlevels(yfac)][yfac]))
p7$scales$scales[[1]]$name <- xlab
p7
```

![plot of chunk colorize-axis-labels](figure/colorize-axis-labels3.png) 



In case you're now wondering what the NMDS/Bray-Curtis plot now looks like on the `openfpp` filtered data.


```r
p5 = plot_ordination(openfpp, ordinate(openfpp, "NMDS", "bray"), color = "bodysite", 
    shape = "GENOTYPE", title = "NMDS on Bray-Curtis distance of filtered open-reference data") + 
    geom_point(size = 5)
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.1478 
## Run 1 stress 0.1448 
## ... New best solution
## ... procrustes: rmse 0.03065  max resid 0.2384 
## Run 2 stress 0.1546 
## Run 3 stress 0.1919 
## Run 4 stress 0.1911 
## Run 5 stress 0.1447 
## ... New best solution
## ... procrustes: rmse 0.001282  max resid 0.006649 
## *** Solution reached
```

```r
p5
```

![plot of chunk NMDS-bray-openfpp](figure/NMDS-bray-openfpp1.png) 

```r
p5 + facet_grid(bodysite ~ GENOTYPE)
```

![plot of chunk NMDS-bray-openfpp](figure/NMDS-bray-openfpp2.png) 



### Network plot

Create network representation of the sample distances. This helps emphasize the two clusters that we see in the UniFrac/MDS plot. One cluster (subnetwork) is comprised entirely of 


```r
openUUFnet = make_network(open, "samples", UniFrac[["unweighted"]][["open"]], 
    max.dist = 0.6, keep.isolates = TRUE)
pnet = plot_network(openUUFnet, open, color = "BODY_PRODUCT", shape = "GENOTYPE", 
    label = NULL)
pnet
```

![plot of chunk plot-network](figure/plot-network1.png) 

```r
pnet + facet_wrap(~GENOTYPE)
```

![plot of chunk plot-network](figure/plot-network2.png) 

```r
plot_network(openUUFnet, open, color = "BODY_SITE")
```

![plot of chunk plot-network](figure/plot-network3.png) 

```r
mtgt = mt(open, get_variable(open, "GENOTYPE"))
```

```
## B=10000
## b=100	b=200	b=300	b=400	b=500	
## b=600	b=700	b=800	b=900	b=1000	b=1100	b=1200	b=1300	b=1400	b=1500	
## b=1600	b=1700	b=1800	b=1900	b=2000	b=2100	b=2200	b=2300	b=2400	b=2500	
## b=2600	b=2700	b=2800	b=2900	b=3000	b=3100	b=3200	b=3300	b=3400	b=3500	
## b=3600	b=3700	b=3800	b=3900	b=4000	b=4100	b=4200	b=4300	b=4400	b=4500	
## b=4600	b=4700	b=4800	b=4900	b=5000	b=5100	b=5200	b=5300	b=5400	b=5500	
## b=5600	b=5700	b=5800	b=5900	b=6000	b=6100	b=6200	b=6300	b=6400	b=6500	
## b=6600	b=6700	b=6800	b=6900	b=7000	b=7100	b=7200	b=7300	b=7400	b=7500	
## b=7600	b=7700	b=7800	b=7900	b=8000	b=8100	b=8200	b=8300	b=8400	b=8500	
## b=8600	b=8700	b=8800	b=8900	b=9000	b=9100	b=9200	b=9300	b=9400	b=9500	
## b=9600	b=9700	b=9800	b=9900	b=10000	r=38	r=76	r=114	r=152	r=190	
## r=228	r=266	r=304	r=342	r=380	r=418	r=456	r=494	r=532	r=570	
## r=608	r=646	r=684	r=722	r=760	r=798	r=836	r=874	r=912	r=950	
## r=988	r=1026	r=1064	r=1102	r=1140	r=1178	r=1216	r=1254	r=1292	r=1330	
## r=1368	r=1406	r=1444	r=1482	r=1520	r=1558	r=1596	r=1634	r=1672	r=1710	
## r=1748	r=1786	r=1824	r=1862	r=1900	r=1938	r=1976	r=2014	r=2052	r=2090	
## r=2128	r=2166	r=2204	r=2242	r=2280	r=2318	r=2356	r=2394	r=2432	r=2470	
## r=2508	r=2546	r=2584	r=2622	r=2660	r=2698	r=2736	r=2774	r=2812	r=2850	
## r=2888	r=2926	r=2964	r=3002	r=3040	r=3078	r=3116	r=3154	r=3192	r=3230	
## r=3268	r=3306	r=3344	r=3382	r=3420	r=3458	r=3496	r=3534	r=3572	r=3610	
## r=3648	r=3686	r=3724	r=3762	r=3800	r=3838	r=3876	
```

```r
mtgt[which(mtgt$adjp < 0.05 & !is.na(mtgt$Phylum)), ]
```

```
##                     index teststat  rawp   adjp plower  Kingdom
## 308076               3852   -4.978 1e-04 0.0182 0.0001 Bacteria
## 212201                276   -4.525 1e-04 0.0182 0.0001 Bacteria
## New.ReferenceOTU61     97    4.261 1e-04 0.0182 0.0001 Bacteria
## 277624                586    4.231 1e-04 0.0182 0.0001 Bacteria
## 463794               2910    3.901 1e-04 0.0182 0.0001 Bacteria
## 175646                397    3.470 1e-04 0.0182 0.0001 Bacteria
## 826624               1472    3.113 1e-04 0.0182 0.0001 Bacteria
## 199118                 57   -4.369 2e-04 0.0365 0.0182 Bacteria
## 374370                319    4.102 2e-04 0.0365 0.0182 Bacteria
## 174056                433    3.897 2e-04 0.0365 0.0182 Bacteria
## 422931               2901    3.835 2e-04 0.0365 0.0182 Bacteria
## 207284                332    3.813 2e-04 0.0365 0.0182 Bacteria
## 278675                 35   -3.771 2e-04 0.0365 0.0182 Bacteria
## 197765               1570    3.767 2e-04 0.0365 0.0182 Bacteria
## 269107               1442    3.635 2e-04 0.0365 0.0182 Bacteria
## 393399               3720    3.532 2e-04 0.0365 0.0182 Bacteria
## New.ReferenceOTU128   229    3.476 2e-04 0.0365 0.0182 Bacteria
## 192222                575    3.336 2e-04 0.0365 0.0182 Bacteria
## 231755               3701    3.327 2e-04 0.0365 0.0182 Bacteria
## New.ReferenceOTU239   895    3.245 2e-04 0.0365 0.0182 Bacteria
## 276149               3608    3.183 2e-04 0.0365 0.0182 Bacteria
## 177352               1363    2.610 2e-04 0.0365 0.0182 Bacteria
##                            Phylum       Class           Order
## 308076              Bacteroidetes Bacteroidia   Bacteroidales
## 212201              Bacteroidetes Bacteroidia   Bacteroidales
## New.ReferenceOTU61  Bacteroidetes Bacteroidia   Bacteroidales
## 277624              Bacteroidetes Bacteroidia   Bacteroidales
## 463794                 Firmicutes     Bacilli Lactobacillales
## 175646              Bacteroidetes Bacteroidia   Bacteroidales
## 826624                 Firmicutes  Clostridia   Clostridiales
## 199118              Bacteroidetes Bacteroidia   Bacteroidales
## 374370              Bacteroidetes Bacteroidia   Bacteroidales
## 174056              Bacteroidetes Bacteroidia   Bacteroidales
## 422931                 Firmicutes     Bacilli Lactobacillales
## 207284              Bacteroidetes Bacteroidia   Bacteroidales
## 278675              Bacteroidetes Bacteroidia   Bacteroidales
## 197765                 Firmicutes  Clostridia   Clostridiales
## 269107                 Firmicutes  Clostridia   Clostridiales
## 393399              Bacteroidetes Bacteroidia   Bacteroidales
## New.ReferenceOTU128 Bacteroidetes Bacteroidia   Bacteroidales
## 192222              Bacteroidetes Bacteroidia   Bacteroidales
## 231755              Bacteroidetes Bacteroidia   Bacteroidales
## New.ReferenceOTU239   Tenericutes  Mollicutes            RF39
## 276149              Bacteroidetes Bacteroidia   Bacteroidales
## 177352                 Firmicutes  Clostridia   Clostridiales
##                                 Family           Genus Species Rank1
## 308076                           S24-7            <NA>    <NA>  <NA>
## 212201                           S24-7            <NA>    <NA>  <NA>
## New.ReferenceOTU61               S24-7            <NA>    <NA>  <NA>
## 277624                  Prevotellaceae      Prevotella    <NA>  <NA>
## 463794                Lactobacillaceae            <NA>    <NA>  <NA>
## 175646                           S24-7            <NA>    <NA>  <NA>
## 826624                 Lachnospiraceae     Coprococcus    <NA>  <NA>
## 199118                           S24-7            <NA>    <NA>  <NA>
## 374370                           S24-7            <NA>    <NA>  <NA>
## 174056                           S24-7            <NA>    <NA>  <NA>
## 422931                Lactobacillaceae            <NA>    <NA>  <NA>
## 207284                           S24-7            <NA>    <NA>  <NA>
## 278675                           S24-7            <NA>    <NA>  <NA>
## 197765                 Lachnospiraceae            <NA>    <NA>  <NA>
## 269107                 Lachnospiraceae  [Ruminococcus]  gnavus  <NA>
## 393399                           S24-7            <NA>    <NA>  <NA>
## New.ReferenceOTU128              S24-7            <NA>    <NA>  <NA>
## 192222                  Prevotellaceae      Prevotella    <NA>  <NA>
## 231755                   Rikenellaceae            <NA>    <NA>  <NA>
## New.ReferenceOTU239               <NA>            <NA>    <NA>  <NA>
## 276149              Porphyromonadaceae Parabacteroides    <NA>  <NA>
## 177352                 Lachnospiraceae            <NA>    <NA>  <NA>
```



It looks like the non-WT genotype never clusters with the WT that was on the left side of the UniFrac/MDS plot shown above. 


---

# Save R image
Save the final R workspace to an `.RData` file for easy re-testing, modification of settings, etc.


```r
save.image("navasetal.RData")
```



---

# References

- Friedrich Leisch. Sweave: Dynamic generation of statistical reports using literate data analysis. In Wolfgang Härdle and Bernd Rönz, editors, Compstat 2002 - Proceedings in Computational Statistics, pages 575-580. Physica Verlag, Heidelberg, 2002. ISBN 3-7908-1517-9

- Yihui Xie (2013) Dynamic Documents with R and knitr. Chapman and Hall/CRC. ISBN 978-1482203530

- Yihui Xie (2013) knitr: A Comprehensive Tool for Reproducible Research in R. In Victoria Stodden, Friedrich Leisch and Roger D. Peng, editors, Implementing Reproducible Computational Research. Chapman and Hall/CRC. ISBN 978-1466561595

- JJ Allaire, Jeffrey Horner, Vicent Marti and Natacha Porte (2013). markdown: Markdown
  rendering for R. R package version 0.5.4. http://CRAN.R-project.org/package=markdown

- Robert Gentleman, et al. Bioconductor: open software development for computational biology and bioinformatics. Genome Biology, 2004. PMID: 15461798

- McMurdie PJ, Holmes S (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE 8(4): e61217.
PMID: 23630581  doi:10.1371/journal.pone.0061217

- Hadley Wickham. ggplot2: Elegant Graphics for Data Analysis. v 6991, Use R! 
Springer, 2009 ISBN 0387981403, 9780387981406
 
- R Core Team. 2013 R: A Language and Environment for Statistical. Vienna, Austria. http://www.R-project.org
```
@Manual{,
  title   = {R: A Language and Environment for Statistical
   Computing},
  author  = {{R Core Team}},
  organization = {R Foundation for Statistical Computing},
  address = {Vienna, Austria},
  year    = 2013,
  url= {http://www.R-project.org}
}
```






