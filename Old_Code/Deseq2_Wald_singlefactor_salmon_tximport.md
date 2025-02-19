```R
#install.packages('BiocManager')
#install.packages('pheatmap')
#install.packages('RColorBrewer')
#install.packages('tidyverse')
#install.packages("devtools") #Also need to install RTools4.0
#install.packages("RColorBrewer")

#library(devtools)
#install_github("cran/lasso2") #Removed from cran

#BiocManager::install("tximportData")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("DEGreport")
#BiocManager::install("clusterProfiler")


#BiocManager::install("apeglm")
#BiocManager::install("EnhancedVolcano")
```


```R
library(tximportData)
library(tximport)
library(DESeq2)
library(tidyverse)
library(DEGreport)
library(pheatmap)
library(RColorBrewer)
library(apeglm)
library(EnhancedVolcano)
```

    Loading required package: S4Vectors
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    
    Attaching package: 'BiocGenerics'
    
    
    The following objects are masked from 'package:stats':
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from 'package:base':
    
        anyDuplicated, append, as.data.frame, basename, cbind, colnames,
        dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
        grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
        union, unique, unsplit, which.max, which.min
    
    
    
    Attaching package: 'S4Vectors'
    
    
    The following objects are masked from 'package:base':
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    
    Attaching package: 'IRanges'
    
    
    The following object is masked from 'package:grDevices':
    
        windows
    
    
    Loading required package: GenomicRanges
    
    Loading required package: GenomeInfoDb
    
    Loading required package: SummarizedExperiment
    
    Loading required package: MatrixGenerics
    
    Loading required package: matrixStats
    
    
    Attaching package: 'MatrixGenerics'
    
    
    The following objects are masked from 'package:matrixStats':
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: 'Biobase'
    
    
    The following object is masked from 'package:MatrixGenerics':
    
        rowMedians
    
    
    The following objects are masked from 'package:matrixStats':
    
        anyMissing, rowMedians
    
    
    -- [1mAttaching core tidyverse packages[22m ------------------------ tidyverse 2.0.0 --
    [32mv[39m [34mdplyr    [39m 1.1.2     [32mv[39m [34mreadr    [39m 2.1.4
    [32mv[39m [34mforcats  [39m 1.0.0     [32mv[39m [34mstringr  [39m 1.5.0
    [32mv[39m [34mggplot2  [39m 3.5.1     [32mv[39m [34mtibble   [39m 3.2.1
    [32mv[39m [34mlubridate[39m 1.9.2     [32mv[39m [34mtidyr    [39m 1.3.0
    [32mv[39m [34mpurrr    [39m 1.0.1     
    -- [1mConflicts[22m ------------------------------------------ tidyverse_conflicts() --
    [31mx[39m [34mlubridate[39m::[32m%within%()[39m masks [34mIRanges[39m::%within%()
    [31mx[39m [34mdplyr[39m::[32mcollapse()[39m     masks [34mIRanges[39m::collapse()
    [31mx[39m [34mdplyr[39m::[32mcombine()[39m      masks [34mBiobase[39m::combine(), [34mBiocGenerics[39m::combine()
    [31mx[39m [34mdplyr[39m::[32mcount()[39m        masks [34mmatrixStats[39m::count()
    [31mx[39m [34mdplyr[39m::[32mdesc()[39m         masks [34mIRanges[39m::desc()
    [31mx[39m [34mtidyr[39m::[32mexpand()[39m       masks [34mS4Vectors[39m::expand()
    [31mx[39m [34mdplyr[39m::[32mfilter()[39m       masks [34mstats[39m::filter()
    [31mx[39m [34mdplyr[39m::[32mfirst()[39m        masks [34mS4Vectors[39m::first()
    [31mx[39m [34mdplyr[39m::[32mlag()[39m          masks [34mstats[39m::lag()
    [31mx[39m [34mggplot2[39m::[32mPosition()[39m   masks [34mBiocGenerics[39m::Position(), [34mbase[39m::Position()
    [31mx[39m [34mpurrr[39m::[32mreduce()[39m       masks [34mGenomicRanges[39m::reduce(), [34mIRanges[39m::reduce()
    [31mx[39m [34mdplyr[39m::[32mrename()[39m       masks [34mS4Vectors[39m::rename()
    [31mx[39m [34mlubridate[39m::[32msecond()[39m   masks [34mS4Vectors[39m::second()
    [31mx[39m [34mlubridate[39m::[32msecond<-()[39m masks [34mS4Vectors[39m::second<-()
    [31mx[39m [34mdplyr[39m::[32mslice()[39m        masks [34mIRanges[39m::slice()
    [36mi[39m Use the conflicted package ([3m[34m<http://conflicted.r-lib.org/>[39m[23m) to force all conflicts to become errors
    Loading required package: ggrepel
    
    Registered S3 methods overwritten by 'ggalt':
      method                  from   
      grid.draw.absoluteGrob  ggplot2
      grobHeight.absoluteGrob ggplot2
      grobWidth.absoluteGrob  ggplot2
      grobX.absoluteGrob      ggplot2
      grobY.absoluteGrob      ggplot2
    
    


```R
citation('DESeq2')
citation('ggplot2')
```


    
      Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change
      and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
      (2014)
    
    A BibTeX entry for LaTeX users is
    
      @Article{,
        title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
        author = {Michael I. Love and Wolfgang Huber and Simon Anders},
        year = {2014},
        journal = {Genome Biology},
        doi = {10.1186/s13059-014-0550-8},
        volume = {15},
        issue = {12},
        pages = {550},
      }
    



    
    To cite ggplot2 in publications, please use
    
      H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
      Springer-Verlag New York, 2016.
    
    A BibTeX entry for LaTeX users is
    
      @Book{,
        author = {Hadley Wickham},
        title = {ggplot2: Elegant Graphics for Data Analysis},
        publisher = {Springer-Verlag New York},
        year = {2016},
        isbn = {978-3-319-24277-4},
        url = {https://ggplot2.tidyverse.org},
      }
    



```R
citation('DEGreport')
citation('pheatmap')
```


    
    To cite package 'DEGreport' in publications use:
    
      Lorena Pantano (2022). DEGreport: Report of DEG analysis. R package
      version 1.30.3. http://lpantano.github.io/DEGreport/
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {DEGreport: Report of DEG analysis},
        author = {Lorena Pantano},
        year = {2022},
        note = {R package version 1.30.3},
        url = {http://lpantano.github.io/DEGreport/},
      }
    



    
    To cite package 'pheatmap' in publications use:
    
      Raivo Kolde (2019). pheatmap: Pretty Heatmaps. R package version
      1.0.12. https://CRAN.R-project.org/package=pheatmap
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {pheatmap: Pretty Heatmaps},
        author = {Raivo Kolde},
        year = {2019},
        note = {R package version 1.0.12},
        url = {https://CRAN.R-project.org/package=pheatmap},
      }
    
    ATTENTION: This citation information has been auto-generated from the
    package DESCRIPTION file and may need manual editing, see
    'help("citation")'.
    



```R
## This workflow is build to work off denovo-transcriptome assembly with trinity & quantificaion with salmon
## Transcript id's , must still be in the format of id_cluster#_gene#_isoform#
```


```R
##Importing metadata will be a slighly individual process.
##The important things to know is to 
##1. we want all cols to be as a factor and not Intergers/flots
##2. We need to make sure we have a column to build a basic model from, i.e. a column that contains the indivial groups
##3. Col 1 needs to be the same was what you ouput salmon folder name is

samples <- read.table(file = "metadata.txt", header = TRUE) ## First we will import our metadata file
samples <- as.data.frame(unclass(samples),stringsAsFactors=TRUE) ## When importing the table we set all string colums to factors
#samples$Time <- as.factor(samples$Time) ## If you have a colum that is numbers, you must manuall set it to a factor

# For this experiment we had one treatment factor across time, so we want a column that groups treatment and time
#samples$Trt_Time <- paste(samples$Treatment, samples$Time, sep="_") #combines treatment column and time column seperate by _
#samples$Trt_Time <- as.factor(samples$Trt_Time) # Setting this new column as a factor
samples ##Inspect your metadata to make sure it looks correct
```


<table class="dataframe">
<caption>A data.frame: 24 × 2</caption>
<thead>
	<tr><th scope=col>SampleName</th><th scope=col>Treatment</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>b5 </td><td>BTV</td></tr>
	<tr><td>b6 </td><td>BTV</td></tr>
	<tr><td>b7 </td><td>BTV</td></tr>
	<tr><td>b8 </td><td>BTV</td></tr>
	<tr><td>b13</td><td>BTV</td></tr>
	<tr><td>b14</td><td>BTV</td></tr>
	<tr><td>b15</td><td>BTV</td></tr>
	<tr><td>b16</td><td>BTV</td></tr>
	<tr><td>v5 </td><td>VSV</td></tr>
	<tr><td>v6 </td><td>VSV</td></tr>
	<tr><td>v7 </td><td>VSV</td></tr>
	<tr><td>v8 </td><td>VSV</td></tr>
	<tr><td>v13</td><td>VSV</td></tr>
	<tr><td>v14</td><td>VSV</td></tr>
	<tr><td>v15</td><td>VSV</td></tr>
	<tr><td>v16</td><td>VSV</td></tr>
	<tr><td>p5 </td><td>PBS</td></tr>
	<tr><td>p6 </td><td>PBS</td></tr>
	<tr><td>p7 </td><td>PBS</td></tr>
	<tr><td>p8 </td><td>PBS</td></tr>
	<tr><td>p13</td><td>PBS</td></tr>
	<tr><td>p14</td><td>PBS</td></tr>
	<tr><td>p15</td><td>PBS</td></tr>
	<tr><td>p16</td><td>PBS</td></tr>
</tbody>
</table>




```R
############################################
# Need to write a function to do this in R #
############################################

## import the table that converts transcript Id to gene ID
## With trinity transcripts this is easy

# Example

# TXNAME	GENEID
# TRINITY_DN179090_c0_g1_i1	TRINITY_DN179090_c0_g1
# TRINITY_DN179085_c0_g1_i1	TRINITY_DN179085_c0_g1
# TRINITY_DN179048_c0_g1_i1	TRINITY_DN179048_c0_g1

tx2gene <- read.table(file = "tx2gene.tsv", header = TRUE) #reading from file
tx2gene
```


<table class="dataframe">
<caption>A data.frame: 22672 × 3</caption>
<thead>
	<tr><th scope=col>transcript_id</th><th scope=col>gene_id</th><th scope=col>gene_name</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>unassigned_transcript_1    </td><td>cson_001691</td><td>cson_001691</td></tr>
	<tr><td>unassigned_transcript_2    </td><td>cson_012105</td><td>cson_012105</td></tr>
	<tr><td>unassigned_transcript_3    </td><td>cson_012105</td><td>cson_012105</td></tr>
	<tr><td>unassigned_transcript_4    </td><td>cson_012105</td><td>cson_012105</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_012105-R1</td><td>cson_012105</td><td>cson_012105</td></tr>
	<tr><td>unassigned_transcript_5    </td><td>cson_002494</td><td>cson_002494</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_014216-R1</td><td>cson_014216</td><td>cson_014216</td></tr>
	<tr><td>unassigned_transcript_6    </td><td>cson_000081</td><td>cson_000081</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_001751-R1</td><td>cson_001751</td><td>cson_001751</td></tr>
	<tr><td>unassigned_transcript_7    </td><td>cson_013771</td><td>cson_013771</td></tr>
	<tr><td>unassigned_transcript_8    </td><td>cson_001250</td><td>cson_001250</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_001249-R1</td><td>cson_001249</td><td>cson_001249</td></tr>
	<tr><td>unassigned_transcript_9    </td><td>cson_000841</td><td>cson_000841</td></tr>
	<tr><td>unassigned_transcript_10   </td><td>cson_001094</td><td>cson_001094</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_014672-R1</td><td>cson_014672</td><td>cson_014672</td></tr>
	<tr><td>unassigned_transcript_11   </td><td>cson_002923</td><td>cson_002923</td></tr>
	<tr><td>unassigned_transcript_12   </td><td>cson_002178</td><td>cson_002178</td></tr>
	<tr><td>unassigned_transcript_13   </td><td>cson_002178</td><td>cson_002178</td></tr>
	<tr><td>unassigned_transcript_14   </td><td>cson_002340</td><td>cson_002340</td></tr>
	<tr><td>unassigned_transcript_15   </td><td>cson_000839</td><td>cson_000839</td></tr>
	<tr><td>unassigned_transcript_16   </td><td>cson_002188</td><td>cson_002188</td></tr>
	<tr><td>unassigned_transcript_17   </td><td>cson_002143</td><td>cson_002143</td></tr>
	<tr><td>unassigned_transcript_18   </td><td>cson_002143</td><td>cson_002143</td></tr>
	<tr><td>unassigned_transcript_19   </td><td>cson_002143</td><td>cson_002143</td></tr>
	<tr><td>unassigned_transcript_20   </td><td>cson_002143</td><td>cson_002143</td></tr>
	<tr><td>unassigned_transcript_21   </td><td>cson_002143</td><td>cson_002143</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_002124-R1</td><td>cson_002124</td><td>cson_002124</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_009393-R1</td><td>cson_009393</td><td>cson_009393</td></tr>
	<tr><td>unassigned_transcript_22   </td><td>cson_002163</td><td>cson_002163</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_002096-R1</td><td>cson_002096</td><td>cson_002096</td></tr>
	<tr><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013559-R2</td><td>cson_013559</td><td>cson_013559</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013559-R3</td><td>cson_013559</td><td>cson_013559</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013559-R1</td><td>cson_013559</td><td>cson_013559</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013559-R4</td><td>cson_013559</td><td>cson_013559</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013559-R5</td><td>cson_013559</td><td>cson_013559</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013564-R1</td><td>cson_013564</td><td>cson_013564</td></tr>
	<tr><td>unassigned_transcript_3879 </td><td>cson_013577</td><td>cson_013577</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013561-R1</td><td>cson_013561</td><td>cson_013561</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013562-R2</td><td>cson_013562</td><td>cson_013562</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013562-R5</td><td>cson_013562</td><td>cson_013562</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013562-R1</td><td>cson_013562</td><td>cson_013562</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013562-R4</td><td>cson_013562</td><td>cson_013562</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013562-R3</td><td>cson_013562</td><td>cson_013562</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013563-R1</td><td>cson_013563</td><td>cson_013563</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013571-R1</td><td>cson_013571</td><td>cson_013571</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013569-R1</td><td>cson_013569</td><td>cson_013569</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013544-R5</td><td>cson_013544</td><td>cson_013544</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013544-R1</td><td>cson_013544</td><td>cson_013544</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013544-R2</td><td>cson_013544</td><td>cson_013544</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013544-R4</td><td>cson_013544</td><td>cson_013544</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013544-R3</td><td>cson_013544</td><td>cson_013544</td></tr>
	<tr><td>unassigned_transcript_3880 </td><td>cson_013609</td><td>cson_013609</td></tr>
	<tr><td>unassigned_transcript_3881 </td><td>cson_013610</td><td>cson_013610</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013543-R1</td><td>cson_013543</td><td>cson_013543</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013549-R1</td><td>cson_013549</td><td>cson_013549</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013547-R1</td><td>cson_013547</td><td>cson_013547</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013547-R2</td><td>cson_013547</td><td>cson_013547</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013548-R1</td><td>cson_013548</td><td>cson_013548</td></tr>
	<tr><td>gnl|WGS:ZZZZ|cson_013545-R1</td><td>cson_013545</td><td>cson_013545</td></tr>
	<tr><td>unassigned_transcript_3882 </td><td>cson_013654</td><td>cson_013654</td></tr>
</tbody>
</table>




```R
## Importing you salmon data
## In this example we have a folder called salmon_quant filled with the output folders from salmon.
## We pull the name of our folders from samples$SampleName, and the quant files from in the folder is called quant.sf

files <- file.path("Quants", samples$SampleName, "quant.sf") ## creating path to each quant files
all(file.exists(files))
names(files) <- samples$SampleName ## associating each quant file with the sample name
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene) #, txOut = TRUE) ## importing quant data, and converting it to gene level
head(txi.salmon$counts)
```


TRUE


    reading in files with read_tsv
    
    1 
    2 
    3 
    4 
    5 
    6 
    7 
    8 
    9 
    10 
    11 
    12 
    13 
    14 
    15 
    16 
    17 
    18 
    19 
    20 
    21 
    22 
    23 
    24 
    
    
    summarizing abundance
    
    summarizing counts
    
    summarizing length
    
    


<table class="dataframe">
<caption>A matrix: 6 × 24 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>b5</th><th scope=col>b6</th><th scope=col>b7</th><th scope=col>b8</th><th scope=col>b13</th><th scope=col>b14</th><th scope=col>b15</th><th scope=col>b16</th><th scope=col>v5</th><th scope=col>v6</th><th scope=col>...</th><th scope=col>v15</th><th scope=col>v16</th><th scope=col>p5</th><th scope=col>p6</th><th scope=col>p7</th><th scope=col>p8</th><th scope=col>p13</th><th scope=col>p14</th><th scope=col>p15</th><th scope=col>p16</th></tr>
</thead>
<tbody>
	<tr><th scope=row>cson_000001</th><td> 83</td><td> 61</td><td> 36</td><td> 59</td><td> 77</td><td> 88</td><td> 96</td><td>101</td><td> 49</td><td> 34</td><td>...</td><td>  99</td><td>104</td><td>159</td><td> 16</td><td> 27</td><td> 67</td><td>118</td><td> 88</td><td> 93</td><td>119</td></tr>
	<tr><th scope=row>cson_000002</th><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>...</td><td>   0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td><td>  0</td></tr>
	<tr><th scope=row>cson_000003</th><td>  2</td><td>  2</td><td>  1</td><td>  3</td><td>  2</td><td>  2</td><td>  3</td><td>  1</td><td>  0</td><td>  1</td><td>...</td><td>   5</td><td>  1</td><td>  2</td><td>  1</td><td>  1</td><td>  3</td><td>  4</td><td>  2</td><td>  0</td><td>  2</td></tr>
	<tr><th scope=row>cson_000004</th><td> 22</td><td> 26</td><td> 12</td><td> 17</td><td>  8</td><td> 21</td><td> 43</td><td> 14</td><td>  8</td><td>  5</td><td>...</td><td>  60</td><td> 25</td><td> 35</td><td> 11</td><td>  2</td><td> 18</td><td> 37</td><td> 23</td><td> 44</td><td> 47</td></tr>
	<tr><th scope=row>cson_000005</th><td> 49</td><td> 39</td><td> 27</td><td> 42</td><td> 43</td><td> 51</td><td> 28</td><td> 58</td><td> 32</td><td> 20</td><td>...</td><td>  86</td><td> 53</td><td>113</td><td> 17</td><td> 13</td><td> 87</td><td> 92</td><td> 68</td><td> 42</td><td> 59</td></tr>
	<tr><th scope=row>cson_000006</th><td>425</td><td>647</td><td>291</td><td>783</td><td>616</td><td>652</td><td>651</td><td>848</td><td>267</td><td>263</td><td>...</td><td>1170</td><td>919</td><td>891</td><td>204</td><td>311</td><td>667</td><td>816</td><td>822</td><td>840</td><td>761</td></tr>
</tbody>
</table>




```R
## Create DESeq2Dataset object
## For a general model you design can the your group column
dds <- DESeqDataSetFromTximport(txi.salmon, colData = samples, design = ~ Treatment)
```

    using counts and average transcript lengths from tximport
    
    


```R
# making sure min size for each group is 3 & filtering our low cout trasncripts (defualt = 10)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
```


```R
dds
```


    class: DESeqDataSet 
    dim: 9403 24 
    metadata(1): version
    assays(2): counts avgTxLength
    rownames(9403): cson_000001 cson_000004 ... cson_014854 cson_014855
    rowData names(0):
    colnames(24): b5 b6 ... p15 p16
    colData names(2): SampleName Treatment



```R
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
```

    using 'avgTxLength' from assays(dds), correcting for library size
    
    


```R
### Plot PCA 
plotPCA(rld, intgroup="Treatment")
```


    
![png](output_12_0.png)
    



```R
# Getting rlog and generating principle comps
# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
```

    using 'avgTxLength' from assays(dds), correcting for library size
    
    


```R
# Create data frame with metadata and PC values for input to ggplot
# creating a dataframe that comtains our PC's and metadata
df_pca <- cbind(samples, pca$x)
df_pca
```


<table class="dataframe">
<caption>A data.frame: 24 × 26</caption>
<thead>
	<tr><th></th><th scope=col>SampleName</th><th scope=col>Treatment</th><th scope=col>PC1</th><th scope=col>PC2</th><th scope=col>PC3</th><th scope=col>PC4</th><th scope=col>PC5</th><th scope=col>PC6</th><th scope=col>PC7</th><th scope=col>PC8</th><th scope=col>...</th><th scope=col>PC15</th><th scope=col>PC16</th><th scope=col>PC17</th><th scope=col>PC18</th><th scope=col>PC19</th><th scope=col>PC20</th><th scope=col>PC21</th><th scope=col>PC22</th><th scope=col>PC23</th><th scope=col>PC24</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>...</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>b5</th><td>b5 </td><td>BTV</td><td> -0.1631254</td><td> -3.8182100</td><td> 13.3559312</td><td> -0.06920559</td><td>  8.1801277</td><td> -1.15513933</td><td> 6.5693381</td><td> -1.4266442</td><td>...</td><td>-5.6087690</td><td> 1.6874211</td><td> 6.6885747</td><td> 3.5624909</td><td>-5.22855991</td><td>  2.5066602</td><td> 2.6373086</td><td>-1.14626074</td><td>-0.57043248</td><td>3.051726e-14</td></tr>
	<tr><th scope=row>b6</th><td>b6 </td><td>BTV</td><td> 15.7003258</td><td>  2.6271222</td><td> -5.6545083</td><td> -0.99149745</td><td>  5.4453932</td><td> -1.21457900</td><td> 2.3373269</td><td>  0.5999971</td><td>...</td><td>-4.0129838</td><td>-2.2199168</td><td>-1.9992393</td><td> 3.1220758</td><td>-0.83636716</td><td>-12.4973138</td><td>-0.6680697</td><td>-2.77247015</td><td>-0.55106906</td><td>2.995520e-14</td></tr>
	<tr><th scope=row>b7</th><td>b7 </td><td>BTV</td><td>  9.8944241</td><td> 11.8809348</td><td> -0.5810218</td><td>  5.58745702</td><td>  3.4429079</td><td>  4.74759195</td><td> 6.1070003</td><td> -3.2882162</td><td>...</td><td> 9.4620291</td><td>-3.3505619</td><td>-1.4011150</td><td>-1.0835328</td><td> 1.63736053</td><td>  2.9336415</td><td>-4.9669508</td><td>-0.51081283</td><td>-1.56151235</td><td>3.040623e-14</td></tr>
	<tr><th scope=row>b8</th><td>b8 </td><td>BTV</td><td> 12.8850490</td><td> 12.5816043</td><td>  2.4907564</td><td> -7.17477367</td><td>  6.6568734</td><td> -1.90348583</td><td>-3.1606894</td><td> -1.1991849</td><td>...</td><td> 0.1687312</td><td> 2.2893580</td><td>-0.9010955</td><td> 2.2164463</td><td>-1.30081396</td><td>  3.4429924</td><td>-1.4583722</td><td>-0.09624037</td><td>-0.82999258</td><td>3.056583e-14</td></tr>
	<tr><th scope=row>b13</th><td>b13</td><td>BTV</td><td>  9.9416491</td><td> 10.5490723</td><td>  4.3231202</td><td> -2.92774752</td><td> -4.8794551</td><td> -1.08778282</td><td>-8.1318648</td><td> -3.7986843</td><td>...</td><td>-0.2026295</td><td> 6.0454607</td><td> 2.1636148</td><td> 2.0187258</td><td> 6.54660829</td><td>  1.4056823</td><td> 3.7804226</td><td> 3.26317456</td><td>-0.60468845</td><td>3.001766e-14</td></tr>
	<tr><th scope=row>b14</th><td>b14</td><td>BTV</td><td>  9.0208123</td><td>  1.0672022</td><td>  5.4544743</td><td> -2.37902205</td><td>-10.4948011</td><td>  9.63426620</td><td> 0.9613898</td><td>  4.8553151</td><td>...</td><td>-1.5557631</td><td>-2.0312753</td><td>-2.5903605</td><td>-1.0586696</td><td>-7.67173433</td><td> -0.1603320</td><td> 1.1273899</td><td> 7.45615531</td><td> 3.49537081</td><td>2.981643e-14</td></tr>
	<tr><th scope=row>b15</th><td>b15</td><td>BTV</td><td>  9.3125849</td><td>-11.3871366</td><td> -4.9149248</td><td> -9.46312539</td><td>  2.5397445</td><td>  4.39859600</td><td> 1.4534654</td><td>  4.0931235</td><td>...</td><td> 0.1511443</td><td>-3.6897429</td><td> 2.5984920</td><td>-6.0681056</td><td> 6.05331166</td><td>  4.7564766</td><td> 3.1343746</td><td>-2.29702459</td><td> 1.57427138</td><td>3.062828e-14</td></tr>
	<tr><th scope=row>b16</th><td>b16</td><td>BTV</td><td>  2.6949752</td><td> -3.3753951</td><td>  1.0342666</td><td>  2.73843197</td><td> -2.0851187</td><td>  0.32053696</td><td>-1.0624251</td><td>  4.6032004</td><td>...</td><td>-2.6213440</td><td>-1.8563341</td><td>-3.3102760</td><td>-3.2704186</td><td> 6.26627872</td><td> -2.8627048</td><td> 1.3853976</td><td>-1.91745254</td><td>-0.63383846</td><td>3.015643e-14</td></tr>
	<tr><th scope=row>v5</th><td>v5 </td><td>VSV</td><td> -9.9733148</td><td> -2.2033100</td><td>  5.5309062</td><td>-10.81710742</td><td>  3.2559331</td><td>-11.09280844</td><td> 0.5459001</td><td>  3.0487373</td><td>...</td><td> 5.9492382</td><td>-2.8799298</td><td>-0.7671140</td><td>-1.4052992</td><td> 0.54552938</td><td> -3.5197120</td><td> 2.5808234</td><td> 3.18401710</td><td> 2.98791910</td><td>3.003153e-14</td></tr>
	<tr><th scope=row>v6</th><td>v6 </td><td>VSV</td><td> -1.4922207</td><td> -0.9426295</td><td> 16.7592052</td><td>  6.14091333</td><td> -4.1710417</td><td>  0.06213773</td><td>-2.0297094</td><td>  3.4843252</td><td>...</td><td> 1.1540381</td><td> 1.1804203</td><td> 3.5352352</td><td>-3.4687093</td><td> 1.52941149</td><td> -0.5995500</td><td>-2.1538678</td><td>-6.83207490</td><td>-1.09089384</td><td>2.951459e-14</td></tr>
	<tr><th scope=row>v7</th><td>v7 </td><td>VSV</td><td> -1.3237795</td><td>  2.8112471</td><td>  2.7299118</td><td>  5.99018604</td><td> -1.6545751</td><td> -3.82247240</td><td> 7.1275724</td><td> -4.5512194</td><td>...</td><td>-5.5470662</td><td>-5.7864939</td><td>-4.4853428</td><td> 1.3735286</td><td> 3.03859241</td><td>  2.3653836</td><td> 6.4329998</td><td> 2.28206749</td><td>-1.55376181</td><td>3.006623e-14</td></tr>
	<tr><th scope=row>v8</th><td>v8 </td><td>VSV</td><td>  0.1473342</td><td> -1.2782382</td><td> -2.8780701</td><td>  7.90723228</td><td> -2.0762805</td><td> -4.07888942</td><td> 0.5906228</td><td>-10.0920546</td><td>...</td><td>-0.5919239</td><td> 0.1064126</td><td>-2.6740797</td><td>-4.3427234</td><td>-2.05080615</td><td>  1.5448245</td><td> 0.1694652</td><td>-3.31613314</td><td> 8.60255135</td><td>2.949724e-14</td></tr>
	<tr><th scope=row>v13</th><td>v13</td><td>VSV</td><td>  7.7628337</td><td> -7.2069999</td><td> -5.6276353</td><td> -4.99903379</td><td> -1.6747257</td><td> -8.51314361</td><td>-0.9829825</td><td>  2.4441158</td><td>...</td><td> 0.2837790</td><td> 1.0544574</td><td>-2.7914004</td><td>-1.6883115</td><td>-7.27799252</td><td>  5.2414015</td><td>-4.0013856</td><td>-3.09720943</td><td>-1.72085867</td><td>3.012174e-14</td></tr>
	<tr><th scope=row>v14</th><td>v14</td><td>VSV</td><td> -0.4928988</td><td> -0.5445463</td><td> -2.2038930</td><td>  0.75503483</td><td> -6.4472995</td><td> -1.80591785</td><td> 4.2693616</td><td> -0.6034894</td><td>...</td><td> 2.5558677</td><td>10.9989382</td><td> 1.0486679</td><td>-4.5009447</td><td>-1.19821904</td><td> -5.7363846</td><td> 0.9544229</td><td>-0.33829507</td><td>-0.68020109</td><td>3.107931e-14</td></tr>
	<tr><th scope=row>v15</th><td>v15</td><td>VSV</td><td>-12.0998395</td><td>  8.3865873</td><td> -2.2168131</td><td> -4.05192501</td><td> -9.0363225</td><td> -3.94544717</td><td> 1.1206563</td><td>  3.8307239</td><td>...</td><td>-6.4210493</td><td>-0.8640172</td><td>-4.2422944</td><td> 3.9246673</td><td> 2.26077246</td><td>  2.8769138</td><td>-5.6589094</td><td>-2.62833953</td><td>-2.08637577</td><td>3.008704e-14</td></tr>
	<tr><th scope=row>v16</th><td>v16</td><td>VSV</td><td> -2.1655492</td><td>  1.3172044</td><td> -2.8504510</td><td>  5.45009841</td><td> -3.6249758</td><td> -8.32339694</td><td>-4.4152538</td><td> -0.6978207</td><td>...</td><td> 4.2911117</td><td>-6.8580082</td><td> 7.7141225</td><td> 1.3051252</td><td>-0.58122031</td><td> -1.6387873</td><td>-2.1798242</td><td> 4.96362521</td><td>-0.56189894</td><td>3.025705e-14</td></tr>
	<tr><th scope=row>p5</th><td>p5 </td><td>PBS</td><td> -9.2329776</td><td> -7.1628934</td><td>  8.9052569</td><td>  2.60990024</td><td> 10.0589032</td><td>  3.51476100</td><td>-6.5238688</td><td>  0.2183055</td><td>...</td><td> 1.1022074</td><td> 2.3711726</td><td>-9.1072677</td><td> 0.5461847</td><td> 1.34903052</td><td> -0.7752126</td><td>-4.0436495</td><td> 3.51732094</td><td> 1.55716702</td><td>3.130482e-14</td></tr>
	<tr><th scope=row>p6</th><td>p6 </td><td>PBS</td><td>-10.3398636</td><td>  4.5652940</td><td>-11.8348917</td><td>  0.96182124</td><td>  4.5151452</td><td>  0.06554464</td><td> 3.5755648</td><td>  1.7384879</td><td>...</td><td> 1.3628972</td><td> 6.5842891</td><td>-1.0361540</td><td> 0.6968428</td><td>-0.02914713</td><td>  2.7538239</td><td> 4.9785035</td><td> 1.54030685</td><td>-0.06707375</td><td>3.064909e-14</td></tr>
	<tr><th scope=row>p7</th><td>p7 </td><td>PBS</td><td> -0.7267678</td><td>  1.7224752</td><td> -8.5821245</td><td> 10.91020455</td><td>  6.6912958</td><td>  0.28673950</td><td>-6.7400994</td><td>  6.4105112</td><td>...</td><td>-6.8098626</td><td> 0.5426859</td><td> 4.4489649</td><td>-4.2824629</td><td>-1.00042674</td><td>  1.2538549</td><td>-1.6492463</td><td> 3.37673841</td><td>-0.20135500</td><td>2.939662e-14</td></tr>
	<tr><th scope=row>p8</th><td>p8 </td><td>PBS</td><td>-11.6945158</td><td>  8.0098950</td><td> -1.0778630</td><td>  1.03875924</td><td>  1.0907916</td><td>  5.51532408</td><td> 8.0155818</td><td>  6.7511449</td><td>...</td><td> 3.6950754</td><td>-0.2828673</td><td> 2.1279638</td><td> 0.6300406</td><td>-0.06478634</td><td> -0.8750005</td><td>-0.5861490</td><td>-0.35444069</td><td>-0.08160819</td><td>2.987801e-14</td></tr>
	<tr><th scope=row>p13</th><td>p13</td><td>PBS</td><td>-15.1418520</td><td>  7.7981757</td><td> -2.9940765</td><td>-10.51736260</td><td>  0.5792555</td><td>  8.62028975</td><td>-4.8765957</td><td> -9.1216180</td><td>...</td><td>-3.3444222</td><td>-3.4327016</td><td> 3.5725892</td><td>-3.3435815</td><td>-2.03562305</td><td> -2.6641658</td><td>-1.2151360</td><td>-2.64365538</td><td> 0.28441138</td><td>2.989275e-14</td></tr>
	<tr><th scope=row>p14</th><td>p14</td><td>PBS</td><td> -2.6936332</td><td> -4.5921040</td><td> -3.3822863</td><td>  4.31908004</td><td> -2.0403019</td><td>  3.97523207</td><td>-8.7432555</td><td>  1.4020966</td><td>...</td><td> 5.7985829</td><td>-2.5010161</td><td>-2.2993714</td><td> 6.1117589</td><td>-4.59864612</td><td>  0.3025253</td><td> 8.0317614</td><td>-5.89698959</td><td>-2.74460506</td><td>3.062481e-14</td></tr>
	<tr><th scope=row>p15</th><td>p15</td><td>PBS</td><td>  1.2845160</td><td>-13.6797084</td><td> -4.4382592</td><td> -0.24796777</td><td> -2.5597048</td><td>  3.52651553</td><td> 1.7854953</td><td> -0.8095036</td><td>...</td><td> 0.3284727</td><td> 2.1770687</td><td> 3.9719263</td><td>10.5829273</td><td> 4.32972799</td><td>  0.4951627</td><td>-4.3800961</td><td>-0.04588152</td><td> 5.80553220</td><td>3.082257e-14</td></tr>
	<tr><th scope=row>p16</th><td>p16</td><td>PBS</td><td> -1.1041663</td><td>-17.1256432</td><td> -1.3470102</td><td> -0.77035092</td><td> -1.7117687</td><td>  2.27552743</td><td> 2.2074686</td><td> -7.8916492</td><td>...</td><td> 0.4126385</td><td> 0.7151804</td><td>-0.2650403</td><td>-1.5780548</td><td> 0.31771932</td><td> -0.5501798</td><td>-2.2512128</td><td> 4.30987462</td><td>-8.76705776</td><td>2.954581e-14</td></tr>
</tbody>
</table>




```R
# Creating a dataframe containing StdDev, proportion explained var and cumulative explained var
eigs <- pca$sdev^2
summary_pca <- rbind(
  SD = sqrt(eigs),
  Proportion = eigs/sum(eigs),
  Cumulative = cumsum(eigs)/sum(eigs))
summary_pca_df <- as.data.frame(summary_pca)
```


```R
# View pca summary
summary_pca_df
```


<table class="dataframe">
<caption>A data.frame: 3 × 24</caption>
<thead>
	<tr><th></th><th scope=col>V1</th><th scope=col>V2</th><th scope=col>V3</th><th scope=col>V4</th><th scope=col>V5</th><th scope=col>V6</th><th scope=col>V7</th><th scope=col>V8</th><th scope=col>V9</th><th scope=col>V10</th><th scope=col>...</th><th scope=col>V15</th><th scope=col>V16</th><th scope=col>V17</th><th scope=col>V18</th><th scope=col>V19</th><th scope=col>V20</th><th scope=col>V21</th><th scope=col>V22</th><th scope=col>V23</th><th scope=col>V24</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>...</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>SD</th><td>8.5174212</td><td>7.8452384</td><td>6.60509309</td><td>5.81433027</td><td>5.32920368</td><td>5.14102673</td><td>4.82304719</td><td>4.63416673</td><td>4.53258125</td><td>4.50551152</td><td>...</td><td>4.11130810</td><td>3.99818287</td><td>3.89869822</td><td>3.83773567</td><td>3.80461265</td><td>3.7652699</td><td>3.62571389</td><td>3.52229845</td><td>3.19308933</td><td>3.064546e-14</td></tr>
	<tr><th scope=row>Proportion</th><td>0.1315484</td><td>0.1116045</td><td>0.07910923</td><td>0.06130115</td><td>0.05149842</td><td>0.04792576</td><td>0.04218055</td><td>0.03894149</td><td>0.03725293</td><td>0.03680929</td><td>...</td><td>0.03064992</td><td>0.02898642</td><td>0.02756186</td><td>0.02670665</td><td>0.02624763</td><td>0.0257076</td><td>0.02383726</td><td>0.02249684</td><td>0.01848806</td><td>1.702949e-30</td></tr>
	<tr><th scope=row>Cumulative</th><td>0.1315484</td><td>0.2431528</td><td>0.32226207</td><td>0.38356322</td><td>0.43506164</td><td>0.48298740</td><td>0.52516795</td><td>0.56410944</td><td>0.60136237</td><td>0.63817167</td><td>...</td><td>0.79996768</td><td>0.82895410</td><td>0.85651596</td><td>0.88322261</td><td>0.90947024</td><td>0.9351778</td><td>0.95901510</td><td>0.98151194</td><td>1.00000000</td><td>1.000000e+00</td></tr>
</tbody>
</table>




```R
# Function to display decimal as percent rounded to 2 sig dig
percent <- function(x, digits = 2, format = "f", ...) {      
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}
```


```R
# Testing percent function
percent(summary_pca_df$V1[2])
```


'13.15%'



```R
####################
# Editable section #
####################

# Set column of data frame to use for shape
shape_col <- 'Treatment'

## Shape list. You need one shape for every level of shape_colum. In the example we have 2 levels
## Link to shapes:
## https://blog.albertkuo.me/post/2019-03-24-point-shapes-in-ggplot_files/figure-html/unnamed-chunk-2-1.png
##  This code will only work well with "BOTH" shapes i.e. 21-25
Shape_list <- list(21, 25)

# Set column of data frame to use for fill color
fill_col <- 'Treatment'

# Set dot size- May need to lower size if you have lots of samples
dot_size <- 7

# plot_dpi
dpi <- 300

# width in Inches of final plot
width <- 6
# heigth in Inches of final plot
heigth <- 6

##########################
# Kinda editable section #
##########################


# The following code draws PCA's for the combination of your first 6 PC's (change at you own risk). 
# I do recommend chaning the title for each one. Ex. ggtitle("PC1 vs. PC2") --> ggtitle("changed title")

pca1 <- ggplot(df_pca) + 
        geom_point(aes(x=PC1, y=PC2, fill = !!ensym(fill_col)), shape = 21, colour="black", size=dot_size, alpha = 8/10) +
        #scale_shape_manual(values=c(Shape_list), labels=c('VSNJV-infected', 'Uninfected')) +
        ggtitle("PC1 vs. PC2") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(face = "bold")) +
        theme(axis.text.x=element_text(size=15)) +
        theme(axis.text.y=element_text(size=15)) +
        guides(fill = guide_legend(override.aes = list(shape=22))) +
        theme(axis.title=element_text(size=14,face="bold")) +
        xlab(paste("PC1 (", percent(summary_pca_df$V1[2]), ")")) +
        ylab(paste("PC2 (", percent(summary_pca_df$V2[2]), ")")) +
        scale_fill_manual(values = c("BTV" = "#ff2c79", "VSV" = "#009093", "PBS" = "#521b92")) +
        theme_bw() +
        labs(shape="Treatment", fill="Treatment") # extra color #FFB347
pca2 <- ggplot(df_pca) + 
        geom_point(aes(x=PC2, y=PC3, fill = !!ensym(fill_col)), shape = 21, colour="black", size=dot_size, alpha = 8/10) +
        #scale_shape_manual(values=c(Shape_list), labels=c('VSNJV-infected', 'Uninfected')) +
        ggtitle("PC2 vs. PC3") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(face = "bold")) +
        theme(axis.text.x=element_text(size=15)) +
        theme(axis.text.y=element_text(size=15)) +
        guides(fill = guide_legend(override.aes = list(shape=22))) +
        theme(axis.title=element_text(size=14,face="bold")) +
        xlab(paste("PC2 (", percent(summary_pca_df$V2[2]), ")")) +
        ylab(paste("PC3 (", percent(summary_pca_df$V3[2]), ")")) +
        scale_fill_manual(values = c("BTV" = "#ff2c79", "VSV" = "#009093", "PBS" = "#521b92")) +
        theme_bw() +
        labs(shape="Infection status", fill="Time (HPI)")
pca3 <- ggplot(df_pca) + 
        geom_point(aes(x=PC3, y=PC4, fill = !!ensym(fill_col)), shape = 21, colour="black", size=dot_size, alpha = 8/10) +
        #scale_shape_manual(values=c(Shape_list), labels=c('VSNJV-infected', 'Uninfected')) +
        ggtitle("PC3 vs. PC4") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(face = "bold")) +
        theme(axis.text.x=element_text(size=15)) +
        theme(axis.text.y=element_text(size=15)) +
        guides(fill = guide_legend(override.aes = list(shape=22))) +
        theme(axis.title=element_text(size=14,face="bold")) +
        xlab(paste("PC3 (", percent(summary_pca_df$V3[2]), ")")) +
        ylab(paste("PC4 (", percent(summary_pca_df$V4[2]), ")")) +
        scale_fill_manual(values = c("BTV" = "#ff2c79", "VSV" = "#009093", "PBS" = "#521b92")) +
        theme_bw() +
        labs(shape="Infection status", fill="Time (HPI)")
pca4 <- ggplot(df_pca) + 
        geom_point(aes(x=PC4, y=PC5, fill = !!ensym(fill_col)), shape = 21, colour="black", size=dot_size, alpha = 8/10) +
        #scale_shape_manual(values=c(Shape_list), labels=c('VSNJV-infected', 'Uninfected')) +
        ggtitle("PC4 vs. PC5") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(face = "bold")) +
        theme(axis.text.x=element_text(size=15)) +
        theme(axis.text.y=element_text(size=15)) +
        guides(fill = guide_legend(override.aes = list(shape=22))) +
        theme(axis.title=element_text(size=14,face="bold")) +
        xlab(paste("PC4 (", percent(summary_pca_df$V4[2]), ")")) +
        ylab(paste("PC5 (", percent(summary_pca_df$V5[2]), ")")) +
        scale_fill_manual(values = c("BTV" = "#ff2c79", "VSV" = "#009093", "PBS" = "#521b92")) +
        theme_bw() +
        labs(shape="Infection status", fill="Time (HPI)")
pca5 <- ggplot(df_pca) + 
        geom_point(aes(x=PC5, y=PC6, fill = !!ensym(fill_col)), shape = 21, colour="black", size=dot_size, alpha = 8/10) +
        #scale_shape_manual(values=c(Shape_list), labels=c('VSNJV-infected', 'Uninfected')) +
        ggtitle("PC5 vs. PC6") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(face = "bold")) +
        theme(axis.text.x=element_text(size=15)) +
        theme(axis.text.y=element_text(size=15)) +
        guides(fill = guide_legend(override.aes = list(shape=22))) +
        theme(axis.title=element_text(size=14,face="bold")) +
        xlab(paste("PC5 (", percent(summary_pca_df$V5[2]), ")")) +
        ylab(paste("PC6 (", percent(summary_pca_df$V6[2]), ")")) +
        scale_fill_manual(values = c("BTV" = "#ff2c79", "VSV" = "#009093", "PBS" = "#521b92")) +
        theme_bw() +
        labs(shape="Infection status", fill="Time (HPI)")
#pca1_3 <- ggplot(df_pca) + 
#        geom_point(aes(x=PC1, y=PC3, fill = !!ensym(fill_col), shape = !!ensym(shape_col)), colour="black", size=dot_size, alpha = 8/10) +
#        scale_shape_manual(values=c(Shape_list), labels=c('VSNJV-infected', 'Uninfected')) +
#        ggtitle("PC1 vs. PC3") +
#        theme(plot.title = element_text(hjust = 0.5)) +
#        theme(plot.title = element_text(face = "bold")) +
#        theme(axis.text.x=element_text(size=15)) +
#        theme(axis.text.y=element_text(size=15)) +
#        guides(fill = guide_legend(override.aes = list(shape=22))) +
#        theme(axis.title=element_text(size=14,face="bold")) +
#        xlab(paste("PC1 (", percent(summary_pca_df$V1[2]), ")")) +
#        ylab(paste("PC3 (", percent(summary_pca_df$V3[2]), ")"))  +
#        scale_fill_manual(values = c("1" = "#ff2c79", "8" = "#009093", "24" = "#521b92", "96" = "#ba82ff")) +
#        theme_bw() +
#        labs(shape="Infection status", fill="Time (HPI)")


###########################
# No reason to Edit below #
###########################


## Making output folder structure

if (!dir.exists("plots")){
  dir.create("plots")
}else{
  print("dir exists")
}

if (!dir.exists("plots/pcas")){
  dir.create("plots/pcas")
}else{
  print("dir exists")
}

## Drawing & Saving the plots

plot(pca1)
png("plots/pcas/PCA_1_2.png",width=width,height=heigth,units="in",res=dpi)
print(pca1)
dev.off()

plot(pca2)
png("plots/pcas/PCA_2_3.png",width=width,height=heigth,units="in",res=dpi)
print(pca2)
dev.off()

plot(pca3)
png("plots/pcas/PCA_3_4.png",width=width,height=heigth,units="in",res=dpi)
print(pca3)
dev.off()

plot(pca4)
png("plots/pcas/PCA_4_5.png",width=width,height=heigth,units="in",res=dpi)
print(pca4)
dev.off()

plot(pca5)
png("plots/pcas/PCA_5_6.png",width=width,height=heigth,units="in",res=dpi)
print(pca5)
dev.off()

#plot(pca1_3)
#png("plots/pcas/PCA_1_3.png",width=width,height=heigth,units="in",res=dpi)
#print(pca1_3)
#dev.off()

## Ouputing PCA matrix & stats
#write.csv(df_pca, "plots/pcas/pca_matrix.csv")
#write.csv(summary_pca_df, "plots/pcas/pca_matrix_summary.csv")
```

    [1] "dir exists"
    [1] "dir exists"
    


<strong>png:</strong> 2



    
![png](output_19_2.png)
    



<strong>png:</strong> 2



    
![png](output_19_4.png)
    



<strong>png:</strong> 2



    
![png](output_19_6.png)
    



<strong>png:</strong> 2



    
![png](output_19_8.png)
    



<strong>png:</strong> 2



    
![png](output_19_10.png)
    



```R
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
```


```R
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## corelation

head(rld_cor)   ## check the output of cor
```


<table class="dataframe">
<caption>A matrix: 6 × 24 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>b5</th><th scope=col>b6</th><th scope=col>b7</th><th scope=col>b8</th><th scope=col>b13</th><th scope=col>b14</th><th scope=col>b15</th><th scope=col>b16</th><th scope=col>v5</th><th scope=col>v6</th><th scope=col>...</th><th scope=col>v15</th><th scope=col>v16</th><th scope=col>p5</th><th scope=col>p6</th><th scope=col>p7</th><th scope=col>p8</th><th scope=col>p13</th><th scope=col>p14</th><th scope=col>p15</th><th scope=col>p16</th></tr>
</thead>
<tbody>
	<tr><th scope=row>b5</th><td>1.0000000</td><td>0.9870490</td><td>0.9873022</td><td>0.9864932</td><td>0.9873576</td><td>0.9885091</td><td>0.9875877</td><td>0.9901022</td><td>0.9893062</td><td>0.9901343</td><td>...</td><td>0.9863463</td><td>0.9885187</td><td>0.9905331</td><td>0.9854352</td><td>0.9868579</td><td>0.9872943</td><td>0.9850676</td><td>0.9885491</td><td>0.9885353</td><td>0.9887959</td></tr>
	<tr><th scope=row>b6</th><td>0.9870490</td><td>1.0000000</td><td>0.9902553</td><td>0.9902149</td><td>0.9899040</td><td>0.9895946</td><td>0.9899625</td><td>0.9891141</td><td>0.9842930</td><td>0.9851973</td><td>...</td><td>0.9840808</td><td>0.9885596</td><td>0.9841408</td><td>0.9859897</td><td>0.9885268</td><td>0.9845208</td><td>0.9819688</td><td>0.9879364</td><td>0.9880681</td><td>0.9864657</td></tr>
	<tr><th scope=row>b7</th><td>0.9873022</td><td>0.9902553</td><td>1.0000000</td><td>0.9896888</td><td>0.9901894</td><td>0.9894028</td><td>0.9858654</td><td>0.9887406</td><td>0.9841648</td><td>0.9862940</td><td>...</td><td>0.9854765</td><td>0.9883084</td><td>0.9841768</td><td>0.9862580</td><td>0.9884731</td><td>0.9868001</td><td>0.9840508</td><td>0.9879253</td><td>0.9858473</td><td>0.9836780</td></tr>
	<tr><th scope=row>b8</th><td>0.9864932</td><td>0.9902149</td><td>0.9896888</td><td>1.0000000</td><td>0.9904680</td><td>0.9885588</td><td>0.9864437</td><td>0.9875349</td><td>0.9849308</td><td>0.9854011</td><td>...</td><td>0.9846567</td><td>0.9872105</td><td>0.9833699</td><td>0.9841043</td><td>0.9861971</td><td>0.9847074</td><td>0.9835593</td><td>0.9860890</td><td>0.9842140</td><td>0.9824907</td></tr>
	<tr><th scope=row>b13</th><td>0.9873576</td><td>0.9899040</td><td>0.9901894</td><td>0.9904680</td><td>1.0000000</td><td>0.9903519</td><td>0.9863307</td><td>0.9891030</td><td>0.9857001</td><td>0.9879888</td><td>...</td><td>0.9868817</td><td>0.9889326</td><td>0.9849518</td><td>0.9848212</td><td>0.9870611</td><td>0.9861545</td><td>0.9853408</td><td>0.9881297</td><td>0.9859122</td><td>0.9843859</td></tr>
	<tr><th scope=row>b14</th><td>0.9885091</td><td>0.9895946</td><td>0.9894028</td><td>0.9885588</td><td>0.9903519</td><td>1.0000000</td><td>0.9893196</td><td>0.9906979</td><td>0.9856492</td><td>0.9898368</td><td>...</td><td>0.9871127</td><td>0.9885491</td><td>0.9860771</td><td>0.9850454</td><td>0.9868470</td><td>0.9866340</td><td>0.9848064</td><td>0.9898039</td><td>0.9889895</td><td>0.9880419</td></tr>
</tbody>
</table>




```R
help(pheatmap)
```



<table width="100%" summary="page for pheatmap {ComplexHeatmap}"><tr><td>pheatmap {ComplexHeatmap}</td><td style="text-align: right;">R Documentation</td></tr></table>

<h2>
Translate pheatmap::pheatmap to ComplexHeatmap::Heatmap
</h2>

<h3>Description</h3>

<p>Translate pheatmap::pheatmap to ComplexHeatmap::Heatmap
</p>


<h3>Usage</h3>

<pre>
pheatmap(mat,
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    kmeans_k = NA,
    breaks = NA,
    border_color = ifelse(nrow(mat) &lt; 100 &amp; ncol(mat) &lt; 100, "grey60", NA),
    cellwidth = NA,
    cellheight = NA,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    clustering_callback = NA,
    cutree_rows = NA,
    cutree_cols = NA,
    treeheight_row = ifelse(class(cluster_rows) == "hclust" || cluster_rows, 50, 0),
    treeheight_col = ifelse(class(cluster_cols) == "hclust" || cluster_cols, 50, 0),
    legend = TRUE,
    legend_breaks = NA,
    legend_labels = NA,
    annotation_row = NA,
    annotation_col = NA,
    annotation = NA,
    annotation_colors = NA,
    annotation_legend = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    drop_levels = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = NA,
    fontsize = 10,
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    angle_col = c("270", "0", "45", "90", "315"),
    display_numbers = FALSE,
    number_format = "%.2f",
    number_color = "grey30",
    fontsize_number = 0.8 * fontsize,
    gaps_row = NULL,
    gaps_col = NULL,
    labels_row = NULL,
    labels_col = NULL,
    filename = NA,
    width = NA,
    height = NA,
    silent = FALSE,
    na_col = "#DDDDDD",
    name = NULL,

    # argument specific for Heatmap()
    heatmap_legend_param = list(),
    ...,
    run_draw = FALSE)
</pre>


<h3>Arguments</h3>

<table summary="R argblock">
<tr valign="top"><td><code>mat</code></td>
<td>
<p>The input matrix.</p>
</td></tr>
<tr valign="top"><td><code>color</code></td>
<td>
<p>The same as in <code>pheatmap</code>. Here you don't necessarily need to generate a long color vector. The discrete colors sent to <code>colorRampPalette</code> are also OK here. E.g. <code>colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)</code> can be simply  replaced as <code>rev(brewer.pal(n = 7, name = "RdYlBu"))</code>.</p>
</td></tr>
<tr valign="top"><td><code>kmeans_k</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>breaks</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>border_color</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>cellwidth</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>cellheight</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>scale</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>cluster_rows</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>cluster_cols</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>clustering_distance_rows</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>clustering_distance_cols</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>clustering_method</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>clustering_callback</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>cutree_rows</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>cutree_cols</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>treeheight_row</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>treeheight_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>legend</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>legend_breaks</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>legend_labels</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation_row</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation_colors</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation_legend</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation_names_row</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>annotation_names_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>drop_levels</code></td>
<td>
<p>Enforced to be <code>TRUE</code>.</p>
</td></tr>
<tr valign="top"><td><code>show_rownames</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>show_colnames</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>main</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>fontsize</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>fontsize_row</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>fontsize_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>angle_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>display_numbers</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>number_format</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>number_color</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>fontsize_number</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>gaps_row</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>gaps_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>labels_row</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>labels_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>filename</code></td>
<td>
<p>Not supported.</p>
</td></tr>
<tr valign="top"><td><code>width</code></td>
<td>
<p>Not supported.</p>
</td></tr>
<tr valign="top"><td><code>height</code></td>
<td>
<p>Not supported.</p>
</td></tr>
<tr valign="top"><td><code>silent</code></td>
<td>
<p>Not supported.</p>
</td></tr>
<tr valign="top"><td><code>na_col</code></td>
<td>
<p>The same as in <code>pheatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>name</code></td>
<td>
<p>Name of the heatmap. This argument is passed to <code>Heatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>heatmap_legend_param</code></td>
<td>
<p>Pass to <code>Heatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>...</code></td>
<td>
<p>Other arguments passed to <code>Heatmap</code>.</p>
</td></tr>
<tr valign="top"><td><code>run_draw</code></td>
<td>
<p>Whether to run <code>draw()</code> function to the heatmap object.</p>
</td></tr>
</table>


<h3>Details</h3>

<p>This function aims to execute <code>pheatmap::pheatmap</code> code purely with ComplexHeatmap.
</p>


<h3>Value</h3>

<p>A <code>Heatmap-class</code> object.
</p>


<h3>See Also</h3>

<p>See <a href="https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/">https://jokergoo.github.io/2020/05/06/translate-from-pheatmap-to-complexheatmap/</a>
</p>
<p><code>compare_pheatmap</code> that compares heatmaps between <code>pheatmap::pheatmap()</code> and <code>ComplexHeatmap::pheatmap()</code>.
</p>


<h3>Examples</h3>

<pre>
# There is no example
NULL

</pre>

<hr /><div style="text-align: center;">[Package <em>ComplexHeatmap</em> version 2.10.0 ]</div>



```R
###Make meta df
meta <- as.data.frame(colData(dds)[,c("Treatment")])
rownames(meta) <- colnames(dds)

colnames(meta) <- c('Treatment')

#meta$`Infection status` <- ifelse(meta$`Infection status` == 'Infected', 'VSNJV-infected', 'Uninfected')

meta

############################
# Change below at own risk #
############################


## make output folder structure

if (!dir.exists("plots/pretty_heatmap")){
  dir.create("plots/pretty_heatmap")
}else{
  print("dir exists")
}

##Function to save heatmap
save_pheatmap_png <- function(x, filename, width=1920, height=1920, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

my_colour = list(
    'Treatment' = c("BTV" = "#ff2c79", "VSV" = "#009093", "PBS" = "#521b92")
)

### Plot heatmap & save

Names_Heatmap <- pheatmap(rld_cor, annotation_colors = my_colour, annotation = meta, show_rownames = TRUE, show_colnames = TRUE, treeheight_row = 0, fontsize = 7)
save_pheatmap_png(Names_Heatmap, "plots/pretty_heatmap/pretty_heatmap_with_samplenames.png")

NoNames_Heatmap <- pheatmap(rld_cor, annotation_colors = my_colour, annotation = meta, show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 0, fontsize = 7)
save_pheatmap_png(NoNames_Heatmap, "plots/pretty_heatmap/pretty_heatmap_without_samplenames.png")
```


<table class="dataframe">
<caption>A data.frame: 24 × 1</caption>
<thead>
	<tr><th></th><th scope=col>Treatment</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>b5</th><td>BTV</td></tr>
	<tr><th scope=row>b6</th><td>BTV</td></tr>
	<tr><th scope=row>b7</th><td>BTV</td></tr>
	<tr><th scope=row>b8</th><td>BTV</td></tr>
	<tr><th scope=row>b13</th><td>BTV</td></tr>
	<tr><th scope=row>b14</th><td>BTV</td></tr>
	<tr><th scope=row>b15</th><td>BTV</td></tr>
	<tr><th scope=row>b16</th><td>BTV</td></tr>
	<tr><th scope=row>v5</th><td>VSV</td></tr>
	<tr><th scope=row>v6</th><td>VSV</td></tr>
	<tr><th scope=row>v7</th><td>VSV</td></tr>
	<tr><th scope=row>v8</th><td>VSV</td></tr>
	<tr><th scope=row>v13</th><td>VSV</td></tr>
	<tr><th scope=row>v14</th><td>VSV</td></tr>
	<tr><th scope=row>v15</th><td>VSV</td></tr>
	<tr><th scope=row>v16</th><td>VSV</td></tr>
	<tr><th scope=row>p5</th><td>PBS</td></tr>
	<tr><th scope=row>p6</th><td>PBS</td></tr>
	<tr><th scope=row>p7</th><td>PBS</td></tr>
	<tr><th scope=row>p8</th><td>PBS</td></tr>
	<tr><th scope=row>p13</th><td>PBS</td></tr>
	<tr><th scope=row>p14</th><td>PBS</td></tr>
	<tr><th scope=row>p15</th><td>PBS</td></tr>
	<tr><th scope=row>p16</th><td>PBS</td></tr>
</tbody>
</table>



    [1] "dir exists"
    


<strong>png:</strong> 2



    
![png](output_23_3.png)
    



<strong>png:</strong> 2



    
![png](output_23_5.png)
    



```R
dds <- DESeq(dds)  ## Running DESeq
```

    estimating size factors
    
    using 'avgTxLength' from assays(dds), correcting for library size
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    
    -- replacing outliers and refitting for 24 genes
    -- DESeq argument 'minReplicatesForReplace' = 7 
    -- original counts are preserved in counts(dds)
    
    estimating dispersions
    
    fitting model and testing
    
    


```R
## Total number of raw counts per sample
colSums(counts(dds))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>b5</dt><dd>2370940</dd><dt>b6</dt><dd>3730238</dd><dt>b7</dt><dd>1910848</dd><dt>b8</dt><dd>3312596</dd><dt>b13</dt><dd>3223154</dd><dt>b14</dt><dd>3750035</dd><dt>b15</dt><dd>4077764</dd><dt>b16</dt><dd>3661860</dd><dt>v5</dt><dd>1545188</dd><dt>v6</dt><dd>1665855</dd><dt>v7</dt><dd>1614800</dd><dt>v8</dt><dd>5353831</dd><dt>v13</dt><dd>2586528</dd><dt>v14</dt><dd>3581171</dd><dt>v15</dt><dd>6338528</dd><dt>v16</dt><dd>4272858</dd><dt>p5</dt><dd>5537743</dd><dt>p6</dt><dd>1022583</dd><dt>p7</dt><dd>1298120</dd><dt>p8</dt><dd>4618357</dd><dt>p13</dt><dd>4992759</dd><dt>p14</dt><dd>4502762</dd><dt>p15</dt><dd>4282605</dd><dt>p16</dt><dd>4498586</dd></dl>




```R
## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>b5</dt><dd>3108547.3629147</dd><dt>b6</dt><dd>3315489.45283171</dd><dt>b7</dt><dd>3200545.2755938</dd><dt>b8</dt><dd>3239652.32003786</dd><dt>b13</dt><dd>3266550.59388752</dd><dt>b14</dt><dd>3098542.80286671</dd><dt>b15</dt><dd>3085760.24053675</dd><dt>b16</dt><dd>3156309.24675634</dd><dt>v5</dt><dd>2859253.29610481</dd><dt>v6</dt><dd>3016716.38937544</dd><dt>v7</dt><dd>3058241.07098637</dd><dt>v8</dt><dd>3203394.84342352</dd><dt>v13</dt><dd>3103867.88215487</dd><dt>v14</dt><dd>3154719.98911154</dd><dt>v15</dt><dd>3216614.7763756</dd><dt>v16</dt><dd>3224924.8033148</dd><dt>p5</dt><dd>2913632.09767826</dd><dt>p6</dt><dd>3053077.01786692</dd><dt>p7</dt><dd>3506367.9266639</dd><dt>p8</dt><dd>3417854.60874683</dd><dt>p13</dt><dd>2861622.26645157</dd><dt>p14</dt><dd>3204363.89353853</dd><dt>p15</dt><dd>3159165.15830354</dd><dt>p16</dt><dd>2955592.0520288</dd></dl>




```R
## Cange figure settings at own risk

## make output folder structure

if (!dir.exists("plots/dispersion_plot")){
  dir.create("plots/dispersion_plot")
}else{
  print("dir exists")
}

##Function to save dispersion plot
save_displot_png <- function(x, filename, width=1600, height=1400, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

## Plot dispersion estimates
plotDispEsts(dds)

## Saving dist plot
save_displot_png(plotDispEsts(dds), "plots/dispersion_plot/dispersion_plot.png")
```

    [1] "dir exists"
    


<strong>png:</strong> 2



    
![png](output_27_2.png)
    



```R
####################################################
# Needs to be rewritten in loops for each contrast #
####################################################

# From here down we start working with indiviual contrasts from the model. Repeat each step for each of your contrasts
# In this example we have 4 contrasts, 1 for each timepoint Infected vs. Uninfected
```


```R
resultsNames(dds)
#Here we can see deseq2 did not run the contrasts we are interested in, so we will specify them anyway
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Intercept'</li><li>'Treatment_PBS_vs_BTV'</li><li>'Treatment_VSV_vs_BTV'</li></ol>




```R
## Define contrasts

contrast_1 <- c("Treatment", "BTV", "PBS")
contrast_8 <- c("Treatment", "VSV", "PBS")

#Extract results table with a given alpha & Log2-Foldchange threshold 
#signifigant results with have a +/-log2FC & adj P-value > set values
#Remeber we are talking in terms of Log2(FC), so 0.58 = 1.5 Fold Change

## Setting values
padj.cutoff <- 0.05 ## Some studies use as high as 0.1 for adjusted p-value, with 0.05 being very common. Can be adjusted lower #sig genes very high
lfc.cutoff <- 0.58 ## Set to 0 if no lfc threshold is desired. To calc actual fold change (FC = 2^lfc)

res_table_C_WC_58FC <- results(dds, contrast=contrast_1, alpha = padj.cutoff, lfcThreshold = lfc.cutoff)
res_table_WC_NP_58FC <- results(dds, contrast=contrast_8, alpha = padj.cutoff, lfcThreshold = lfc.cutoff)
```


```R
#plotting significantly differentially expressed genes (in blue)
plotMA(res_table_C_WC_58FC, ylim=c(-3,3), main='res_table_C_WC_1FC')
plotMA(res_table_WC_NP_58FC, ylim=c(-3,3), main='res_table_WC_NP_1FC')
```


    
![png](output_31_0.png)
    



    
![png](output_31_1.png)
    



```R
print('summary of res_table_C_WC_58FC')
summary(res_table_C_WC_58FC, alpha = 0.05)
print('summary of res_table_WC_NP_58FC')
summary(res_table_WC_NP_58FC, alpha = 0.05)
```

    [1] "summary of res_table_C_WC_58FC"
    
    out of 9402 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0.58 (up)    : 3, 0.032%
    LFC < -0.58 (down) : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.011%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    [1] "summary of res_table_WC_NP_58FC"
    
    out of 9402 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0.58 (up)    : 0, 0%
    LFC < -0.58 (down) : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.011%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    


```R
## Creating tables of only signficiant results

# table 1-hour
res_table_C_WC_58FC_tb <- res_table_C_WC_58FC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_res_C_WC_58FC_pval <- filter(res_table_C_WC_58FC_tb, padj < padj.cutoff)

# table 8-hour
res_table_WC_NP_58FC_tb <- res_table_WC_NP_58FC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_res_WC_NP_58FC_pval <- filter(res_table_WC_NP_58FC_tb, padj < padj.cutoff)
```


```R
print("Number Significant reults abouve 2FC in C Vs. WC")
nrow(sig_res_C_WC_58FC_pval)
print("Number Significant reults abouve 2FC in WC Vs. NP")
nrow(sig_res_WC_NP_58FC_pval)
```

    [1] "Number Significant reults abouve 2FC in C Vs. WC"
    


3


    [1] "Number Significant reults abouve 2FC in WC Vs. NP"
    


0



```R
## Define contrasts

contrast_1 <- c("Treatment", "BTV", "PBS")
contrast_8 <- c("Treatment", "VSV", "PBS")

#Extract results table with a given alpha & Log2-Foldchange threshold 
#signifigant results with have a +/-log2FC & adj P-value > set values
#Remeber we are talking in terms of Log2(FC), so 0.58 = 1.5 Fold Change

## Setting values
padj.cutoff <- 0.05 ## Some studies use as high as 0.1 for adjusted p-value, with 0.05 being very common. Can be adjusted lower #sig genes very high
lfc.cutoff <- 0 ## Set to 0 if no lfc threshold is desired. To calc actual fold change (FC = 2^lfc)

res_table_C_WC_0FC <- results(dds, contrast=contrast_1, alpha = padj.cutoff, lfcThreshold = lfc.cutoff)
res_table_WC_NP_0FC <- results(dds, contrast=contrast_8, alpha = padj.cutoff, lfcThreshold = lfc.cutoff)
```


```R
#plotting significantly differentially expressed genes (in blue)
plotMA(res_table_C_WC_0FC, ylim=c(-3,3), main='res_table_C_WC_0FC')
plotMA(res_table_WC_NP_0FC, ylim=c(-3,3), main='res_table_WC_NP_0FC')
```


    
![png](output_36_0.png)
    



    
![png](output_36_1.png)
    



```R
print('summary of res_table_C_WC_0FC')
summary(res_table_C_WC_0FC, alpha = 0.05)
print('summary of res_table_WC_NP_0FC')
summary(res_table_WC_NP_0FC, alpha = 0.05)
```

    [1] "summary of res_table_C_WC_0FC"
    
    out of 9402 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 224, 2.4%
    LFC < 0 (down)     : 231, 2.5%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.011%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    [1] "summary of res_table_WC_NP_0FC"
    
    out of 9402 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0 (up)       : 25, 0.27%
    LFC < 0 (down)     : 5, 0.053%
    outliers [1]       : 0, 0%
    low counts [2]     : 1, 0.011%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    


```R
## Creating tables of only signficiant results

# table 1-hour
res_table_C_WC_0FC_tb <- res_table_C_WC_0FC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_res_C_WC_0FC_pval <- filter(res_table_C_WC_0FC_tb, padj < padj.cutoff)

# table 8-hour
res_table_WC_NP_0FC_tb <- res_table_WC_NP_0FC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_res_WC_NP_0FC_pval <- filter(res_table_WC_NP_0FC_tb, padj < padj.cutoff)
```


```R
print("Number Significant reults abouve 0FC in C Vs. WC")
nrow(sig_res_C_WC_0FC_pval)
print("Number Significant reults abouve 0FC in WC Vs. NP")
nrow(sig_res_WC_NP_0FC_pval)
```

    [1] "Number Significant reults abouve 0FC in C Vs. WC"
    


455


    [1] "Number Significant reults abouve 0FC in WC Vs. NP"
    


30



```R
##Outputting data

## make output folder structure

if (!dir.exists("deseq2_results")){
  dir.create("deseq2_results")
}else{
  print("dir exists")
}

if (!dir.exists("deseq2_results/sig_results")){
  dir.create("deseq2_results/sig_results_58FC")
}else{
  print("dir exists")
}

if (!dir.exists("deseq2_results/sig_results")){
  dir.create("deseq2_results/sig_results_0FC")
}else{
  print("dir exists")
}

if (!dir.exists("deseq2_results/results")){
  dir.create("deseq2_results/results")
}else{
  print("dir exists")
}

if (!dir.exists("deseq2_results/norm_counts")){
  dir.create("deseq2_results/norm_counts")
}else{
  print("dir exists")
}

if (!dir.exists("deseq2_results/cluster_profile")){
  dir.create("deseq2_results/cluster_profile")
}else{
  print("dir exists")
}

##Significant results (0.05 & 0.58lfc)
write.csv(sig_res_C_WC_0FC_pval, "deseq2_results/sig_results_0FC/sig_BTV_PBS_0FC.csv")
write.csv(sig_res_WC_NP_0FC_pval, "deseq2_results/sig_results_0FC/sig_VSV_PBS_0FC.csv")

##Significant results (0.05 & 0.58lfc)
write.csv(sig_res_C_WC_58FC_pval, "deseq2_results/sig_results_58FC/sig_BTV_PBS_1FC.csv")
write.csv(sig_res_WC_NP_58FC_pval, "deseq2_results/sig_results_58FC/sig_VSV_PBS_1FC.csv")
#write.csv(sig_res_WC_M_58FC_pval, "deseq2_results/sig_results_58FC/sig_WC_M_1FC.csv")
#write.csv(sig_res_NP_M_58FC_pval, "deseq2_results/sig_results_58FC/sig_NP_M_1FC.csv")

## Full results
#write.csv(res_table_C_WC_58FC, "deseq2_results/results/result_C_WC.csv")
#write.csv(res_table_WC_NP_58FC, "deseq2_results/results/result_WC_NP.csv")
#write.csv(res_table_WC_M_58FC, "deseq2_results/results/result_WC_M.csv")
#write.csv(res_table_NP_M_58FC, "deseq2_results/results/result_NP_M.csv")

## Normalized Count Data - avgTxLength (median of rations- ONlY FOR INTER-SAMPLE GENE COMPARSION, NOT SUTABLE FOR INTRA-SAMPLES GENE COMPARISION) 
write.csv(assays(dds)[["avgTxLength"]], "deseq2_results/norm_counts/Normalized_Counts_avgTxLength.csv")
```

    [1] "dir exists"
    

    Warning message in dir.create("deseq2_results/sig_results_58FC"):
    "'deseq2_results\sig_results_58FC' already exists"
    Warning message in dir.create("deseq2_results/sig_results_0FC"):
    "'deseq2_results\sig_results_0FC' already exists"
    

    [1] "dir exists"
    [1] "dir exists"
    [1] "dir exists"
    

    Warning message in file(file, ifelse(append, "a", "w")):
    "cannot open file 'deseq2_results/sig_results_0FC/sig_BTV_PBS_0FC.csv': Permission denied"
    


    Error in file(file, ifelse(append, "a", "w")): cannot open the connection
    Traceback:
    

    1. write.csv(sig_res_C_WC_0FC_pval, "deseq2_results/sig_results_0FC/sig_BTV_PBS_0FC.csv")

    2. eval.parent(Call)

    3. eval(expr, p)

    4. eval(expr, p)

    5. utils::write.table(sig_res_C_WC_0FC_pval, "deseq2_results/sig_results_0FC/sig_BTV_PBS_0FC.csv", 
     .     col.names = NA, sep = ",", dec = ".", qmethod = "double")

    6. file(file, ifelse(append, "a", "w"))



```R
assays(dds)
```


    List of length 8
    names(8): counts avgTxLength ... replaceCounts replaceCooks



```R

```
