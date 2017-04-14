DESeq2 output README
====================

This document describes:
1) the DESeq2 output structure
2) output files and plots



Example output
=======
This document gives and example directory/file structure.
```
$ tree deseq2_diffex/
deseq2_diffex/
├── counts
│   ├── depth_normalized_counts.txt
│   ├── raw_counts.txt
│   └── rlog_normalized_counts.txt
├── gene_lists
│   ├── diet
│   │   ├── HF_v_DRG.txt
│   │   ├── HF_v_ND.txt
│   │   └── ND_v_DRG.txt
... ... ...
│   └── male.diet
│       └── HF_v_ND.txt
└── plots
    ├── summary_plots
    │   ├── BoxPlot_DepthNormalizedCounts.html
    │   ├── BoxPlot_RawCounts.html
    │   ├── BoxPlot_RlogNormalizedCounts.html
    │   ├── BoxPlot.pdf
    │   ├── CorrelMatrix.pdf
    │   ├── CorrelMatrix_All.html
    │   ├── CorrelMatrix_Fe^HF.html
    │   ├── CorrelMatrix_Me^HF.html
    │   ├── DensityPlot.pdf
    │   ├── DensityPlot_DepthNormalizedCounts.html
    │   ├── DensityPlot_RawCounts.html
    │   ├── DensityPlot_RlogNormalizedCounts.html
    │   ├── Dispersion.pdf
    │   ├── PCAplot_All.html
    │   ├── PCAplot_All.pdf
    │   ├── Heatmap_Samples.pdf
    │   ├── Heatmap_TopExp.pdf
    │   └── Heatmap_TopVar.pdf
    └── comparison_plots
        ├── diet
        │   ├── MAplot_HF_v_DRG.html
        │   ├── MAplot_HF_v_DRG.pdf
        │   ├── MAplot_HF_v_ND.html
        │   ├── MAplot_HF_v_ND.pdf
        │   ├── MAplot_ND_v_DRG.html
        │   ├── MAplot_ND_v_DRG.pdf
        │   ├── PCAplot.html
        │   ├── PCAplot.pdf
        │   ├── VolcanoPlot_HF_v_DRG.html
        │   ├── VolcanoPlot_HF_v_DRG.pdf
        │   ├── VolcanoPlot_HF_v_ND.html
        │   ├── VolcanoPlot_HF_v_ND.pdf
        │   ├── VolcanoPlot_ND_v_DRG.html
        │   └── VolcanoPlot_ND_v_DRG.pdf
        ...
        └── male.diet
            ├── MAplot_HF_v_ND.html
            ├── MAplot_HF_v_ND.pdf
            ├── PCAplot.html
            ├── PCAplot.pdf
            ├── Volcano_HF_v_ND.html
            └── Volcano_HF_v_ND.pdf
```


General structure
====================
The output directory (*deseq2_diffex* in this example) contains 3 sub-directories:  
**`deseq2_diffex/`**  
- **`counts/`** - contains count data  
- **`gene_lists/`** - contains differential expression data  
- **`plots/`** - contains plots  


`counts/` directory structure
===============
**`counts/`**  
- **`depth_normalized_counts.txt`** - raw counts normalized by read depth per sample  
- **`raw_counts.txt`** - raw counts per sample, used in DESeq2 differential expression calculations  
- **`rlog_normalized_counts.txt`** - rlog normalized counts (log2), used for plots in spotting batch effects, etc.  

_The_ *_counts.txt are tab-delimited, and contain the a single 'id' column, with the gene id, and 'Sample_XXX' for each sample._

```
$ cat deseq2_diffex/counts/depth_normalized_counts.txt | head | awk 'BEGIN {FS = OFS = "\t"} {print $1, $2, $3}'
id	Sample_61483	Sample_61484
4930547N16Rik	33.3268545915621	29.6157351920312
6330406I15Rik	233.287982140935	359.876054909228
6330527O06Rik	1107.360486656	1983.35681134512
AU021092	97.7082782343526	98.7191173067708
Adora2a	273.431693353498	70.0008286357102
Aldh1a3	48.475424860454	17.9489304194129
Alox12	65.1388521562351	63.7187029889157
Alox12b	250.708837950161	230.643755889455
Alox12e	4.54457108066756	3.58978608388257
```

`gene_lists directory/` structure
==================
**`gene_lists/`**
- **`phenotype_directory/`** - a phenotype of interest (i.e., 'diet', 'cell', 'gender')  
   - **`condition_v_control.txt`** - differential expression data file for the comparison 'condition' vs 'control' (i.e., samples where phenotype 'gender' is 'male' vs samples where phenotype 'gender' is 'female')

_The differential expression *.txt files are tab-delimited and contain:_

```
$ head deseq2_diffex/gene_lists/diet/HF_v_DRG.txt
id	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	Condition	Control	Call
Egr2	458.291	2.894	1.025	2.838	0.004	0.360	HF	DRG NO
Fos	2545.962	2.813	1.026	2.757	0.005	0.360	HF	DRG NO
```
Column descriptions:  
- `id` - gene identifier
- `baseMean` - mean of normalized counts for all samples
- `log2FoldChange` - log2 fold change value for condition vs control
- `lfcSE` - standard error of the log2 fold change for condition vs control
- `stat` - Wald statistic for condition vs control
- `pvalue` - p-value
- `padj` - Benjamini-Hochberg adjusted p-value
- `Condition` - status value considered 'condition', used as numerator in fold change calculations 
- `Control` - status value considered 'control', used as denominator in fold change calculations
- `Call` - whether gene is considered differentially expressed (YES or NO) based upon input log2FoldChange and padj cutoffs (must pass both cutoffs)

plots directory structure
=============

The plots in the plots directory fall into 2 sub-directories:

**`plots/`**
- **`summary_plots/`** - plots **not specific** to any factor and comparison, generally summarizing all samples    
	- `Box plots` - displays expression data (counts) stats with a boxplot, useful for spotting sample specific issues (raw, depth-normalized, and rlog-normalized counts)
  	- `Correlation matrix plots` - spearman correlation plot matrices for all samples, or samples of a specific group (using rlog-normalized counts)
    - `Density plots` - histogram plots showing the count data distribution for each sample, useful for spotting sample specific issues (raw, depth-normalized, and rlog-normalized counts)
    - `Dispersion plots` - shows the gene-wise relative variability of true expression (between biological replicates) estimates (black), the fitted values (red), and the final maximum a posteriori estimates used in testing (blue)
    - `PCA plots` - plot of the first two principal components (along which the variation in gene expression is maximal for all samples) to view the variation within between samples, all samples are labeled
    - `Heatmaps` - heatmap and dendrogram showing the
    	- Euclidean distance between samples (Heatmap_Samples)
    	- 500 most highly expressed genes (Heatmap_TopExp)
    	- 500 most higly variably expressed genes (Heatmap_TopVar)  
- **`comparison_plots/`** - plots **specific** to a phenotype and comparison
    - **`phenotype_directory/`** - corresponding to each phenotype of interest (**diet**, **male.diet**)  
        - `MA plots` - scatter plot that visualizes, for a specific comparison, the relationship between mean normalized expression  (x-axis) and log2 fold-change (y-axis) values for each gene. Dots represent genes, differentially regulated genes (red, blue respectively) _pass both_ the log2 fold-change and adjusted PValue cutoff specified in the config file. Non-significantly regulated genes (gray) _fail at least one_ of these filters.
        - `PCA plots` - plot of the first two principal components (along which the variation in gene expression is maximal for all samples) to view the variation within between samples, only samples within the specific comparison are labeled
        - `Volcano plots` - scatter plot that visualizes, for a specific comparison, the relationship between log2 fold-change (x-axis) and significance (y-axis) values for each gene. Dots represent genes, differentially regulated genes (red, blue respectively) _pass both_ the log2 fold-change and adjusted PValue cutoff specified in the config file. Non-significantly regulated genes (gray) _fail at least one_ of these filters.

