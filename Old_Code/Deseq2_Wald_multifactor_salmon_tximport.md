```R
#BiocManager::install("apeglm")
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
    [32mv[39m [34mggplot2  [39m 3.4.2     [32mv[39m [34mtibble   [39m 3.2.1
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
    


```R
samples <- read.table(file = "metadata_numerictime.tsv.txt", header = TRUE)
samples <- as.data.frame(unclass(samples),stringsAsFactors=TRUE)
samples$Time <- as.factor(samples$Time)
samples$Trt_Time <- paste(samples$Treatment, samples$Time, sep="_")
samples$Trt_Time <- as.factor(samples$Trt_Time)
samples
```


<table class="dataframe">
<caption>A data.frame: 48 × 4</caption>
<thead>
	<tr><th scope=col>SampleName</th><th scope=col>Treatment</th><th scope=col>Time</th><th scope=col>Trt_Time</th></tr>
	<tr><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>C_1_A </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><td>C_1_B </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><td>C_1_C </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><td>C_1_D </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><td>C_1_E </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><td>C_1_F </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><td>C_8_A </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><td>C_8_B </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><td>C_8_C </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><td>C_8_D </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><td>C_8_E </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><td>C_8_F </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><td>C_24_A</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><td>C_24_B</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><td>C_24_C</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><td>C_24_D</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><td>C_24_E</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><td>C_24_F</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><td>C_96_A</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><td>C_96_B</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><td>C_96_C</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><td>C_96_D</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><td>C_96_E</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><td>C_96_F</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><td>V_1_A </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><td>V_1_B </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><td>V_1_C </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><td>V_1_D </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><td>V_1_E </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><td>V_1_F </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><td>V_8_A </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><td>V_8_B </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><td>V_8_C </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><td>V_8_D </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><td>V_8_E </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><td>V_8_F </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><td>V_24_A</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><td>V_24_B</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><td>V_24_C</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><td>V_24_D</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><td>V_24_E</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><td>V_24_F</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><td>V_96_A</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><td>V_96_B</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><td>V_96_C</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><td>V_96_D</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><td>V_96_E</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><td>V_96_F</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
</tbody>
</table>




```R
txdb <- read.table(file = "gene_2_tx.tsv", header = TRUE)
tx2gene <- txdb
```


```R
files <- file.path("salmon_quant", "salmon", samples$SampleName, "quant.sf")
names(files) <- samples$SampleName
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)
```

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
    25 
    26 
    27 
    28 
    29 
    30 
    31 
    32 
    33 
    34 
    35 
    36 
    37 
    38 
    39 
    40 
    41 
    42 
    43 
    44 
    45 
    46 
    47 
    48 
    
    
    summarizing abundance
    
    summarizing counts
    
    summarizing length
    
    


<table class="dataframe">
<caption>A matrix: 6 × 48 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>C_1_A</th><th scope=col>C_1_B</th><th scope=col>C_1_C</th><th scope=col>C_1_D</th><th scope=col>C_1_E</th><th scope=col>C_1_F</th><th scope=col>C_8_A</th><th scope=col>C_8_B</th><th scope=col>C_8_C</th><th scope=col>C_8_D</th><th scope=col>...</th><th scope=col>V_24_C</th><th scope=col>V_24_D</th><th scope=col>V_24_E</th><th scope=col>V_24_F</th><th scope=col>V_96_A</th><th scope=col>V_96_B</th><th scope=col>V_96_C</th><th scope=col>V_96_D</th><th scope=col>V_96_E</th><th scope=col>V_96_F</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TRINITY_DN0_c0_g1</th><td>4490.708</td><td>3898.111</td><td>6816.250</td><td>5244.265</td><td>5658.865</td><td>4545.606</td><td>5102.081</td><td>5131.738</td><td>5683.426</td><td>5264.201</td><td>...</td><td>6451.364</td><td>5375.294</td><td>6658.806</td><td>5705.084</td><td>3785.423</td><td>3262.685</td><td>3030.182</td><td>3394.920</td><td>3033.362</td><td>2928.931</td></tr>
	<tr><th scope=row>TRINITY_DN0_c0_g2</th><td>   4.000</td><td>   6.000</td><td>  17.000</td><td>   8.000</td><td>  10.000</td><td>   6.000</td><td>   8.000</td><td>  22.000</td><td>  12.000</td><td>  11.000</td><td>...</td><td>   6.000</td><td>  11.000</td><td>  12.000</td><td>   6.000</td><td>   3.000</td><td>   6.000</td><td>   3.000</td><td>   4.000</td><td>   8.000</td><td>   4.000</td></tr>
	<tr><th scope=row>TRINITY_DN0_c0_g3</th><td>   0.000</td><td>   0.000</td><td>   0.000</td><td>   1.000</td><td>   0.000</td><td>   0.000</td><td>   1.000</td><td>   1.000</td><td>   0.000</td><td>   1.000</td><td>...</td><td>   0.000</td><td>   0.000</td><td>   1.000</td><td>   0.000</td><td>   0.000</td><td>   0.000</td><td>   0.000</td><td>   1.000</td><td>   2.000</td><td>   0.000</td></tr>
	<tr><th scope=row>TRINITY_DN0_c1_g1</th><td>  73.250</td><td>  64.943</td><td> 152.050</td><td> 124.758</td><td> 127.939</td><td>  93.405</td><td> 453.622</td><td> 373.047</td><td> 509.331</td><td> 431.117</td><td>...</td><td> 219.163</td><td> 183.695</td><td> 250.006</td><td> 195.472</td><td> 394.395</td><td> 306.893</td><td> 265.258</td><td> 344.766</td><td> 249.255</td><td> 263.474</td></tr>
	<tr><th scope=row>TRINITY_DN0_c10_g1</th><td>1388.850</td><td>1023.287</td><td>1871.033</td><td>1460.069</td><td>1574.582</td><td>1358.728</td><td>1479.769</td><td>1413.415</td><td>1599.582</td><td>1450.988</td><td>...</td><td>1699.625</td><td>1406.224</td><td>1852.171</td><td>1526.316</td><td>1613.952</td><td>1437.908</td><td>1305.970</td><td>1533.899</td><td>1257.485</td><td>1239.962</td></tr>
	<tr><th scope=row>TRINITY_DN0_c11_g1</th><td>   5.311</td><td>   4.000</td><td>   5.000</td><td>   2.000</td><td>   5.000</td><td>   6.000</td><td>   6.000</td><td>   6.000</td><td>   8.170</td><td>  11.000</td><td>...</td><td>   9.000</td><td>   4.000</td><td>   9.000</td><td>   4.002</td><td>   7.000</td><td>   3.000</td><td>   4.098</td><td>   4.000</td><td>   5.000</td><td>   6.000</td></tr>
</tbody>
</table>




```R
## Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi.salmon, colData = samples, design = ~ Time + Treatment)
```

    using counts and average transcript lengths from tximport
    
    


```R
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
```


```R
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
```

    rlog() may take a few minutes with 30 or more samples,
    vst() is a much faster transformation
    
    using 'avgTxLength' from assays(dds), correcting for library size
    
    


```R
### Plot PCA 
plotPCA(rld, intgroup="Treatment")
```


    
![png](output_8_0.png)
    



```R
### Plot PCA 
plotPCA(rld, intgroup="Time")
```


    
![png](output_9_0.png)
    



```R
# Input is a matrix of log transformed values
rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
```

    rlog() may take a few minutes with 30 or more samples,
    vst() is a much faster transformation
    
    using 'avgTxLength' from assays(dds), correcting for library size
    
    


```R
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(samples, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = Time))
```


    
![png](output_11_0.png)
    



```R
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
```


```R
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function

head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
```


<table class="dataframe">
<caption>A matrix: 6 × 48 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>C_1_A</th><th scope=col>C_1_B</th><th scope=col>C_1_C</th><th scope=col>C_1_D</th><th scope=col>C_1_E</th><th scope=col>C_1_F</th><th scope=col>C_8_A</th><th scope=col>C_8_B</th><th scope=col>C_8_C</th><th scope=col>C_8_D</th><th scope=col>...</th><th scope=col>V_24_C</th><th scope=col>V_24_D</th><th scope=col>V_24_E</th><th scope=col>V_24_F</th><th scope=col>V_96_A</th><th scope=col>V_96_B</th><th scope=col>V_96_C</th><th scope=col>V_96_D</th><th scope=col>V_96_E</th><th scope=col>V_96_F</th></tr>
</thead>
<tbody>
	<tr><th scope=row>C_1_A</th><td>1.0000000</td><td>0.9981413</td><td>0.9982691</td><td>0.9981970</td><td>0.9981102</td><td>0.9981103</td><td>0.9954881</td><td>0.9957069</td><td>0.9952711</td><td>0.9954089</td><td>...</td><td>0.9900807</td><td>0.9899364</td><td>0.9905971</td><td>0.9904765</td><td>0.9809202</td><td>0.9799958</td><td>0.9806420</td><td>0.9803470</td><td>0.9811947</td><td>0.9807837</td></tr>
	<tr><th scope=row>C_1_B</th><td>0.9981413</td><td>1.0000000</td><td>0.9983720</td><td>0.9981165</td><td>0.9981235</td><td>0.9980838</td><td>0.9955503</td><td>0.9957781</td><td>0.9954102</td><td>0.9954751</td><td>...</td><td>0.9902365</td><td>0.9900278</td><td>0.9907457</td><td>0.9905742</td><td>0.9808168</td><td>0.9799028</td><td>0.9805062</td><td>0.9803309</td><td>0.9811738</td><td>0.9807059</td></tr>
	<tr><th scope=row>C_1_C</th><td>0.9982691</td><td>0.9983720</td><td>1.0000000</td><td>0.9983522</td><td>0.9983784</td><td>0.9983029</td><td>0.9956466</td><td>0.9957883</td><td>0.9955113</td><td>0.9955289</td><td>...</td><td>0.9902015</td><td>0.9899673</td><td>0.9907470</td><td>0.9905796</td><td>0.9804467</td><td>0.9796102</td><td>0.9802553</td><td>0.9800712</td><td>0.9808014</td><td>0.9803739</td></tr>
	<tr><th scope=row>C_1_D</th><td>0.9981970</td><td>0.9981165</td><td>0.9983522</td><td>1.0000000</td><td>0.9982431</td><td>0.9982720</td><td>0.9955352</td><td>0.9957200</td><td>0.9953022</td><td>0.9954682</td><td>...</td><td>0.9899531</td><td>0.9898123</td><td>0.9905354</td><td>0.9904441</td><td>0.9803518</td><td>0.9794133</td><td>0.9800676</td><td>0.9798040</td><td>0.9806913</td><td>0.9802333</td></tr>
	<tr><th scope=row>C_1_E</th><td>0.9981102</td><td>0.9981235</td><td>0.9983784</td><td>0.9982431</td><td>1.0000000</td><td>0.9982928</td><td>0.9953964</td><td>0.9955862</td><td>0.9951607</td><td>0.9952639</td><td>...</td><td>0.9899586</td><td>0.9897422</td><td>0.9906383</td><td>0.9905126</td><td>0.9804850</td><td>0.9797088</td><td>0.9803151</td><td>0.9800550</td><td>0.9808670</td><td>0.9804887</td></tr>
	<tr><th scope=row>C_1_F</th><td>0.9981103</td><td>0.9980838</td><td>0.9983029</td><td>0.9982720</td><td>0.9982928</td><td>1.0000000</td><td>0.9955438</td><td>0.9957169</td><td>0.9953247</td><td>0.9954869</td><td>...</td><td>0.9903439</td><td>0.9901841</td><td>0.9908985</td><td>0.9908038</td><td>0.9811215</td><td>0.9802781</td><td>0.9809057</td><td>0.9806640</td><td>0.9814852</td><td>0.9811099</td></tr>
</tbody>
</table>




```R
###Make meta df
meta <- as.data.frame(colData(dds)[,c("Treatment","Time")])
### Plot heatmap
pheatmap(rld_cor, annotation = meta)
```


    
![png](output_14_0.png)
    



```R
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
  		fontsize_row = 10, height=20)
```


    
![png](output_15_0.png)
    



```R
dds <- DESeq(dds)
```

    estimating size factors
    
    using 'avgTxLength' from assays(dds), correcting for library size
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    
    


```R
## Check the size factors
##sizeFactors(dds) Not needed as we used tximport
```


```R
## Total number of raw counts per sample
colSums(counts(dds))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>C_1_A</dt><dd>36055406</dd><dt>C_1_B</dt><dd>30446178</dd><dt>C_1_C</dt><dd>53806552</dd><dt>C_1_D</dt><dd>40477853</dd><dt>C_1_E</dt><dd>41892415</dd><dt>C_1_F</dt><dd>35113843</dd><dt>C_8_A</dt><dd>43966924</dd><dt>C_8_B</dt><dd>42490565</dd><dt>C_8_C</dt><dd>50225755</dd><dt>C_8_D</dt><dd>44979520</dd><dt>C_8_E</dt><dd>38386529</dd><dt>C_8_F</dt><dd>44702168</dd><dt>C_24_A</dt><dd>40003875</dd><dt>C_24_B</dt><dd>42831251</dd><dt>C_24_C</dt><dd>39030598</dd><dt>C_24_D</dt><dd>32109636</dd><dt>C_24_E</dt><dd>49968229</dd><dt>C_24_F</dt><dd>39001001</dd><dt>C_96_A</dt><dd>36680925</dd><dt>C_96_B</dt><dd>33524817</dd><dt>C_96_C</dt><dd>38435406</dd><dt>C_96_D</dt><dd>31432991</dd><dt>C_96_E</dt><dd>33741279</dd><dt>C_96_F</dt><dd>37315021</dd><dt>V_1_A</dt><dd>45070930</dd><dt>V_1_B</dt><dd>38081505</dd><dt>V_1_C</dt><dd>36479058</dd><dt>V_1_D</dt><dd>43700443</dd><dt>V_1_E</dt><dd>50061276</dd><dt>V_1_F</dt><dd>38109419</dd><dt>V_8_A</dt><dd>34289698</dd><dt>V_8_B</dt><dd>32478084</dd><dt>V_8_C</dt><dd>38089026</dd><dt>V_8_D</dt><dd>32647532</dd><dt>V_8_E</dt><dd>34524733</dd><dt>V_8_F</dt><dd>36862025</dd><dt>V_24_A</dt><dd>44793808</dd><dt>V_24_B</dt><dd>41951252</dd><dt>V_24_C</dt><dd>41413513</dd><dt>V_24_D</dt><dd>34456690</dd><dt>V_24_E</dt><dd>42774843</dd><dt>V_24_F</dt><dd>35999234</dd><dt>V_96_A</dt><dd>40914371</dd><dt>V_96_B</dt><dd>37440735</dd><dt>V_96_C</dt><dd>32943516</dd><dt>V_96_D</dt><dd>38711099</dd><dt>V_96_E</dt><dd>33549261</dd><dt>V_96_F</dt><dd>31594934</dd></dl>




```R
## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>C_1_A</dt><dd>41113885.3495362</dd><dt>C_1_B</dt><dd>41008970.6067055</dd><dt>C_1_C</dt><dd>40496791.4154581</dd><dt>C_1_D</dt><dd>40142190.3648165</dd><dt>C_1_E</dt><dd>39725311.5185737</dd><dt>C_1_F</dt><dd>40134902.4197954</dd><dt>C_8_A</dt><dd>38943454.5109683</dd><dt>C_8_B</dt><dd>39213477.2763253</dd><dt>C_8_C</dt><dd>39817255.6409423</dd><dt>C_8_D</dt><dd>38989017.0244476</dd><dt>C_8_E</dt><dd>37635946.3132214</dd><dt>C_8_F</dt><dd>37549212.158686</dd><dt>C_24_A</dt><dd>36360746.5353238</dd><dt>C_24_B</dt><dd>36770755.7836036</dd><dt>C_24_C</dt><dd>35192669.6580488</dd><dt>C_24_D</dt><dd>37017226.444471</dd><dt>C_24_E</dt><dd>36325273.2351795</dd><dt>C_24_F</dt><dd>35170713.7121779</dd><dt>C_96_A</dt><dd>39626222.7112004</dd><dt>C_96_B</dt><dd>40935483.5681873</dd><dt>C_96_C</dt><dd>38850225.092282</dd><dt>C_96_D</dt><dd>39571200.5304722</dd><dt>C_96_E</dt><dd>39920121.0132916</dd><dt>C_96_F</dt><dd>39530406.8822178</dd><dt>V_1_A</dt><dd>39779471.0664924</dd><dt>V_1_B</dt><dd>39355745.7795022</dd><dt>V_1_C</dt><dd>39586168.1745887</dd><dt>V_1_D</dt><dd>40061091.920706</dd><dt>V_1_E</dt><dd>39429501.3539472</dd><dt>V_1_F</dt><dd>39815785.6037011</dd><dt>V_8_A</dt><dd>37494385.6264681</dd><dt>V_8_B</dt><dd>38096734.9477661</dd><dt>V_8_C</dt><dd>38138175.8483135</dd><dt>V_8_D</dt><dd>39357054.2243405</dd><dt>V_8_E</dt><dd>38479681.9710815</dd><dt>V_8_F</dt><dd>38170946.5540468</dd><dt>V_24_A</dt><dd>34960975.8630387</dd><dt>V_24_B</dt><dd>36663099.2146325</dd><dt>V_24_C</dt><dd>36583668.2655217</dd><dt>V_24_D</dt><dd>37428299.0063126</dd><dt>V_24_E</dt><dd>36119392.2450068</dd><dt>V_24_F</dt><dd>36152449.9912563</dd><dt>V_96_A</dt><dd>40401938.6339858</dd><dt>V_96_B</dt><dd>41069667.5428085</dd><dt>V_96_C</dt><dd>40301280.8097273</dd><dt>V_96_D</dt><dd>40736910.5592872</dd><dt>V_96_E</dt><dd>40943561.9102875</dd><dt>V_96_F</dt><dd>40784617.9277906</dd></dl>




```R
## Plot dispersion estimates
plotDispEsts(dds)
```


    
![png](output_20_0.png)
    



```R
?results
```



<table width="100%" summary="page for results {DESeq2}"><tr><td>results {DESeq2}</td><td style="text-align: right;">R Documentation</td></tr></table>

<h2>Extract results from a DESeq analysis</h2>

<h3>Description</h3>

<p><code>results</code> extracts a result table from a DESeq analysis giving base means across samples,
log2 fold changes, standard errors, test statistics, p-values and adjusted p-values;
<code>resultsNames</code> returns the names of the estimated effects (coefficents) of the model;
<code>removeResults</code> returns a <code>DESeqDataSet</code> object with results columns removed.
</p>


<h3>Usage</h3>

<pre>
results(
  object,
  contrast,
  name,
  lfcThreshold = 0,
  altHypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
  listValues = c(1, -1),
  cooksCutoff,
  independentFiltering = TRUE,
  alpha = 0.1,
  filter,
  theta,
  pAdjustMethod = "BH",
  filterFun,
  format = c("DataFrame", "GRanges", "GRangesList"),
  saveCols = NULL,
  test,
  addMLE = FALSE,
  tidy = FALSE,
  parallel = FALSE,
  BPPARAM = bpparam(),
  minmu = 0.5
)

resultsNames(object)

removeResults(object)
</pre>


<h3>Arguments</h3>

<table summary="R argblock">
<tr valign="top"><td><code>object</code></td>
<td>
<p>a DESeqDataSet, on which one
of the following functions has already been called:
<code>DESeq</code>, <code>nbinomWaldTest</code>, or <code>nbinomLRT</code></p>
</td></tr>
<tr valign="top"><td><code>contrast</code></td>
<td>
<p>this argument specifies what comparison to extract from
the <code>object</code> to build a results table. one of either:
</p>

<ul>
<li><p> a character vector with exactly three elements:
the name of a factor in the design formula,
the name of the numerator level for the fold change,
and the name of the denominator level for the fold change
(simplest case)
</p>
</li>
<li><p> a list of 2 character vectors: the names of the fold changes
for the numerator, and the names of the fold changes
for the denominator.
these names should be elements of <code>resultsNames(object)</code>.
if the list is length 1, a second element is added which is the
empty character vector, <code>character()</code>.
(more general case, can be to combine interaction terms and main effects)
</p>
</li>
<li><p> a numeric contrast vector with one element
for each element in <code>resultsNames(object)</code> (most general case)
</p>
</li></ul>

<p>If specified, the <code>name</code> argument is ignored.</p>
</td></tr>
<tr valign="top"><td><code>name</code></td>
<td>
<p>the name of the individual effect (coefficient) for
building a results table. Use this argument rather than <code>contrast</code>
for continuous variables, individual effects or for individual interaction terms.
The value provided to <code>name</code> must be an element of <code>resultsNames(object)</code>.</p>
</td></tr>
<tr valign="top"><td><code>lfcThreshold</code></td>
<td>
<p>a non-negative value which specifies a log2 fold change
threshold. The default value is 0, corresponding to a test that
the log2 fold changes are equal to zero. The user can
specify the alternative hypothesis using the <code>altHypothesis</code> argument,
which defaults to testing
for log2 fold changes greater in absolute value than a given threshold.
If <code>lfcThreshold</code> is specified,
the results are for Wald tests, and LRT p-values will be overwritten.</p>
</td></tr>
<tr valign="top"><td><code>altHypothesis</code></td>
<td>
<p>character which specifies the alternative hypothesis,
i.e. those values of log2 fold change which the user is interested in
finding. The complement of this set of values is the null hypothesis which
will be tested. If the log2 fold change specified by <code>name</code>
or by <code>contrast</code> is written as <i> beta </i>, then the possible values for
<code>altHypothesis</code> represent the following alternate hypotheses:
</p>

<ul>
<li><p> greaterAbs: <i> |beta| &gt; lfcThreshold </i>,
and p-values are two-tailed
</p>
</li>
<li><p> lessAbs: <i> |beta| &lt; lfcThreshold </i>,
p-values are the maximum of the upper and lower tests.
The Wald statistic given is positive, an SE-scaled distance from the closest boundary
</p>
</li>
<li><p> greater: <i> beta &gt; lfcThreshold </i>
</p>
</li>
<li><p> less: <i> beta &lt; -lfcThreshold </i>
</p>
</li></ul>
</td></tr>
<tr valign="top"><td><code>listValues</code></td>
<td>
<p>only used if a list is provided to <code>contrast</code>:
a numeric of length two: the log2 fold changes in the list are multiplied by these values.
the first number should be positive and the second negative. 
by default this is <code>c(1,-1)</code></p>
</td></tr>
<tr valign="top"><td><code>cooksCutoff</code></td>
<td>
<p>theshold on Cook's distance, such that if one or more
samples for a row have a distance higher, the p-value for the row is
set to NA. The default cutoff is the .99 quantile of the F(p, m-p) distribution,
where p is the number of coefficients being fitted and m is the number of samples.
Set to <code>Inf</code> or <code>FALSE</code> to disable the resetting of p-values to NA.
Note: this test excludes the Cook's distance of samples belonging to experimental
groups with only 2 samples.</p>
</td></tr>
<tr valign="top"><td><code>independentFiltering</code></td>
<td>
<p>logical, whether independent filtering should be
applied automatically</p>
</td></tr>
<tr valign="top"><td><code>alpha</code></td>
<td>
<p>the significance cutoff used for optimizing the independent
filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a
value other than 0.1, <code>alpha</code> should be set to that value.</p>
</td></tr>
<tr valign="top"><td><code>filter</code></td>
<td>
<p>the vector of filter statistics over which the independent
filtering will be optimized. By default the mean of normalized counts is used.</p>
</td></tr>
<tr valign="top"><td><code>theta</code></td>
<td>
<p>the quantiles at which to assess the number of rejections
from independent filtering</p>
</td></tr>
<tr valign="top"><td><code>pAdjustMethod</code></td>
<td>
<p>the method to use for adjusting p-values, see <code>?p.adjust</code></p>
</td></tr>
<tr valign="top"><td><code>filterFun</code></td>
<td>
<p>an optional custom function for performing independent filtering
and p-value adjustment, with arguments <code>res</code> (a DESeqResults object),
<code>filter</code> (the quantitity for filtering tests),
<code>alpha</code> (the target FDR),
<code>pAdjustMethod</code>. This function should return a DESeqResults object
with a <code>padj</code> column.</p>
</td></tr>
<tr valign="top"><td><code>format</code></td>
<td>
<p>character, either <code>"DataFrame"</code>,
<code>"GRanges"</code>, or <code>"GRangesList"</code>,
whether the results should be printed as a <code>DESeqResults</code> DataFrame,
or if the results DataFrame should be attached as metadata columns to
the <code>GRanges</code> or <code>GRangesList</code> <code>rowRanges</code> of the <code>DESeqDataSet</code>.
If the <code>rowRanges</code> is a <code>GRangesList</code>, and <code>GRanges</code> is requested, 
the range of each gene will be returned</p>
</td></tr>
<tr valign="top"><td><code>saveCols</code></td>
<td>
<p>character or numeric vector, the columns of
<code>mcols(object)</code> to pass into the <code>results</code> output</p>
</td></tr>
<tr valign="top"><td><code>test</code></td>
<td>
<p>this is automatically detected internally if not provided.
the one exception is after <code>nbinomLRT</code> has been run, <code>test="Wald"</code>
will generate Wald statistics and Wald test p-values.</p>
</td></tr>
<tr valign="top"><td><code>addMLE</code></td>
<td>
<p>if <code>betaPrior=TRUE</code> was used (non-default),
this logical argument specifies if the &quot;unshrunken&quot; maximum likelihood estimates (MLE)
of log2 fold change should be added as a column to the results table (default is FALSE).
This argument is preserved for backward compatability, as now <code>betaPrior=TRUE</code>
by default and the recommended pipeline is
to generate shrunken MAP estimates using <code>lfcShrink</code>.
This argument functionality is only implemented for <code>contrast</code>
specified as three element character vectors.</p>
</td></tr>
<tr valign="top"><td><code>tidy</code></td>
<td>
<p>whether to output the results table with rownames as a first column 'row'.
the table will also be coerced to <code>data.frame</code></p>
</td></tr>
<tr valign="top"><td><code>parallel</code></td>
<td>
<p>if FALSE, no parallelization. if TRUE, parallel
execution using <code>BiocParallel</code>, see next argument <code>BPPARAM</code></p>
</td></tr>
<tr valign="top"><td><code>BPPARAM</code></td>
<td>
<p>an optional parameter object passed internally
to <code>bplapply</code> when <code>parallel=TRUE</code>.
If not specified, the parameters last registered with
<code>register</code> will be used.</p>
</td></tr>
<tr valign="top"><td><code>minmu</code></td>
<td>
<p>lower bound on the estimated count (used when calculating contrasts)</p>
</td></tr>
</table>


<h3>Details</h3>

<p>The results table when printed will provide the information about
the comparison, e.g. &quot;log2 fold change (MAP): condition treated vs untreated&quot;, meaning
that the estimates are of log2(treated / untreated), as would be returned by
<code>contrast=c("condition","treated","untreated")</code>.
Multiple results can be returned for analyses beyond a simple two group comparison,
so <code>results</code> takes arguments <code>contrast</code> and <code>name</code> to help
the user pick out the comparisons of interest for printing a results table.
The use of the <code>contrast</code> argument is recommended for exact specification
of the levels which should be compared and their order.
</p>
<p>If <code>results</code> is run without specifying <code>contrast</code> or <code>name</code>,
it will return the comparison of the last level of the last variable in the
design formula over the first level of this variable. For example, for a simple two-group
comparison, this would return the log2 fold changes of the second group over the
first group (the reference level). Please see examples below and in the vignette. 
</p>
<p>The argument <code>contrast</code> can be used to generate results tables for
any comparison of interest, for example, the log2 fold change between
two levels of a factor, and its usage is described below. It can also
accomodate more complicated numeric comparisons.
Note that <code>contrast</code> will set to 0 the estimated LFC in a
comparison of two groups, where all of the counts in the two groups
are equal to 0 (while other groups have positive counts), while
<code>name</code> will not automatically set these LFC to 0.
The test statistic used for a contrast is:
</p>
<p style="text-align: center;"><i> c' beta / sqrt( c' Sigma c ) </i></p>

<p>The argument <code>name</code> can be used to generate results tables for
individual effects, which must be individual elements of <code>resultsNames(object)</code>.
These individual effects could represent continuous covariates, effects
for individual levels, or individual interaction effects.
</p>
<p>Information on the comparison which was used to build the results table,
and the statistical test which was used for p-values (Wald test or likelihood ratio test)
is stored within the object returned by <code>results</code>. This information is in
the metadata columns of the results table, which is accessible by calling <code>mcols</code>
on the <code>DESeqResults</code> object returned by <code>results</code>.
</p>
<p>On p-values:
</p>
<p>By default, independent filtering is performed to select a set of genes
for multiple test correction which maximizes the number of adjusted
p-values less than a given critical value <code>alpha</code> (by default 0.1).
See the reference in this man page for details on independent filtering.
The filter used for maximizing the number of rejections is the mean
of normalized counts for all samples in the dataset.
Several arguments from the <code>filtered_p</code> function of
the genefilter package (used within the <code>results</code> function)
are provided here to control the independent filtering behavior.
In DESeq2 version &gt;= 1.10, the threshold that is chosen is
the lowest quantile of the filter for which the
number of rejections is close to the peak of a curve fit
to the number of rejections over the filter quantiles.
'Close to' is defined as within 1 residual standard deviation.
The adjusted p-values for the genes which do not pass the filter threshold
are set to <code>NA</code>. 
</p>
<p>By default, <code>results</code> assigns a p-value of <code>NA</code>
to genes containing count outliers, as identified using Cook's distance.
See the <code>cooksCutoff</code> argument for control of this behavior.
Cook's distances for each sample are accessible as a matrix &quot;cooks&quot;
stored in the <code>assays()</code> list. This measure is useful for identifying rows where the
observed counts might not fit to a Negative Binomial distribution.
</p>
<p>For analyses using the likelihood ratio test (using <code>nbinomLRT</code>),
the p-values are determined solely by the difference in deviance between
the full and reduced model formula. A single log2 fold change is printed
in the results table for consistency with other results table outputs,
however the test statistic and p-values may nevertheless involve
the testing of one or more log2 fold changes.
Which log2 fold change is printed in the results table can be controlled
using the <code>name</code> argument, or by default this will be the estimated
coefficient for the last element of <code>resultsNames(object)</code>.
</p>
<p>If <code>useT=TRUE</code> was specified when running <code>DESeq</code> or <code>nbinomWaldTest</code>,
then the p-value generated by <code>results</code> will also make use of the
t distribution for the Wald statistic, using the degrees of freedom
in <code>mcols(object)$tDegreesFreedom</code>.
</p>


<h3>Value</h3>

<p>For <code>results</code>: a <code>DESeqResults</code> object, which is
a simple subclass of DataFrame. This object contains the results columns:
<code>baseMean</code>, <code>log2FoldChange</code>, <code>lfcSE</code>, <code>stat</code>,
<code>pvalue</code> and <code>padj</code>,
and also includes metadata columns of variable information.
The <code>lfcSE</code> gives the standard error of the <code>log2FoldChange</code>.
For the Wald test, <code>stat</code> is the Wald statistic: the <code>log2FoldChange</code>
divided by <code>lfcSE</code>, which is compared to a standard Normal distribution
to generate a two-tailed <code>pvalue</code>. For the likelihood ratio test (LRT),
<code>stat</code> is the difference in deviance between the reduced model and the full model,
which is compared to a chi-squared distribution to generate a <code>pvalue</code>.
</p>
<p>For <code>resultsNames</code>: the names of the columns available as results,
usually a combination of the variable name and a level
</p>
<p>For <code>removeResults</code>: the original <code>DESeqDataSet</code> with results metadata columns removed
</p>


<h3>References</h3>

<p>Richard Bourgon, Robert Gentleman, Wolfgang Huber: Independent
filtering increases detection power for high-throughput experiments.
PNAS (2010), <a href="http://dx.doi.org/10.1073/pnas.0914005107">http://dx.doi.org/10.1073/pnas.0914005107</a>
</p>


<h3>See Also</h3>

<p><code>DESeq</code>, <code>lfcShrink</code>
</p>


<h3>Examples</h3>

<pre>

## Example 1: two-group comparison

dds &lt;- makeExampleDESeqDataSet(m=4)

dds &lt;- DESeq(dds)
res &lt;- results(dds, contrast=c("condition","B","A"))

# with more than two groups, the call would look similar, e.g.:
# results(dds, contrast=c("condition","C","A"))
# etc.

## Example 2: two conditions, two genotypes, with an interaction term

dds &lt;- makeExampleDESeqDataSet(n=100,m=12)
dds$genotype &lt;- factor(rep(rep(c("I","II"),each=3),2))

design(dds) &lt;- ~ genotype + condition + genotype:condition
dds &lt;- DESeq(dds) 
resultsNames(dds)

# the condition effect for genotype I (the main effect)
results(dds, contrast=c("condition","B","A"))

# the condition effect for genotype II
# this is, by definition, the main effect *plus* the interaction term
# (the extra condition effect in genotype II compared to genotype I).
results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))

# the interaction term, answering: is the condition effect *different* across genotypes?
results(dds, name="genotypeII.conditionB")

## Example 3: two conditions, three genotypes

# ~~~ Using interaction terms ~~~

dds &lt;- makeExampleDESeqDataSet(n=100,m=18)
dds$genotype &lt;- factor(rep(rep(c("I","II","III"),each=3),2))
design(dds) &lt;- ~ genotype + condition + genotype:condition
dds &lt;- DESeq(dds)
resultsNames(dds)

# the condition effect for genotype I (the main effect)
results(dds, contrast=c("condition","B","A"))

# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
# (the extra condition effect in genotype III compared to genotype I).
results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))

# the interaction term for condition effect in genotype III vs genotype I.
# this tests if the condition effect is different in III compared to I
results(dds, name="genotypeIII.conditionB")

# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))

# Note that a likelihood ratio could be used to test if there are any
# differences in the condition effect between the three genotypes.

# ~~~ Using a grouping variable ~~~

# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.

dds$group &lt;- factor(paste0(dds$genotype, dds$condition))
design(dds) &lt;- ~ group
dds &lt;- DESeq(dds)
resultsNames(dds)

# the condition effect for genotypeIII
results(dds, contrast=c("group", "IIIB", "IIIA"))

</pre>

<hr /><div style="text-align: center;">[Package <em>DESeq2</em> version 1.34.0 ]</div>



```R
resultsNames(dds)

# Not interested in contrasting any of these
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Intercept'</li><li>'Time_8_vs_1'</li><li>'Time_24_vs_1'</li><li>'Time_96_vs_1'</li><li>'Treatment_Uninfected_vs_Infected'</li></ol>




```R
## From here down dosent work
```


```R
## Define contrasts, extract results table, and shrink the log2 fold changes

contrast_1 <- c("Time1.TreatmentInfected", "Time1.TreatmentUninfected")
contrast_8 <- c("Time8.TreatmentInfected", "Time8.TreatmentUninfected")
contrast_24 <- c("Time24.TreatmentInfected", "Time24.TreatmentUninfected")
contrast_96 <- c("Time96.TreatmentInfected", "Time96.TreatmentUninfected")

res_table_1 <- results(dds, contrast=contrast_1, alpha = 0.05, lfcThreshold = 0.58)
res_table_8 <- results(dds, contrast=contrast_8, alpha = 0.05, lfcThreshold = 0.58)
res_table_24 <- results(dds, contrast=contrast_24, alpha = 0.05, lfcThreshold = 0.58)
res_table_96 <- results(dds, contrast=contrast_96, alpha = 0.05, lfcThreshold = 0.58)
```


    Error in checkContrast(contrast, resNames): 'contrast', as a character vector of length 3, should have the form:
    contrast = c('factorName','numeratorLevel','denominatorLevel'),
    see the manual page of ?results for more information
    Traceback:
    

    1. results(dds, contrast = contrast_1, alpha = 0.05, lfcThreshold = 0.58)

    2. checkContrast(contrast, resNames)

    3. stop("'contrast', as a character vector of length 3, should have the form:\ncontrast = c('factorName','numeratorLevel','denominatorLevel'),\nsee the manual page of ?results for more information")



```R
class(res_table_1)
```


'DESeqResults'



```R
mcols(res_table_1, use.names=T)
```


    DataFrame with 6 rows and 2 columns
                           type            description
                    <character>            <character>
    baseMean       intermediate mean of normalized c..
    log2FoldChange      results log2 fold change (ML..
    lfcSE               results standard error: Trt_..
    stat                results Wald statistic: Trt_..
    pvalue              results Wald test p-value: T..
    padj                results   BH adjusted p-values



```R
res_table_1 %>% data.frame() %>% View()
```


<table class="dataframe">
<caption>A data.frame: 17779 × 6</caption>
<thead>
	<tr><th></th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>TRINITY_DN0_c0_g1</th><td> 4898.712104</td><td> 0.210116016</td><td>0.02331269</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c0_g2</th><td>    9.299546</td><td> 0.438524391</td><td>0.34326004</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c1_g1</th><td>  259.542125</td><td>-0.367983662</td><td>0.10516348</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c10_g1</th><td> 1446.017687</td><td>-0.166869963</td><td>0.03026031</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c11_g1</th><td>    6.576095</td><td> 0.753797609</td><td>0.41938029</td><td> 0.4144153</td><td>0.67856998</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c142_g1</th><td>   84.534552</td><td> 0.196676960</td><td>0.12715418</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c15_g1</th><td>   71.669703</td><td> 0.011562284</td><td>0.14955774</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c2_g1</th><td>    5.545240</td><td>-1.530496463</td><td>0.48734815</td><td>-1.9503439</td><td>0.05113515</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c20_g1</th><td>   35.304042</td><td> 0.537921960</td><td>0.20142297</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c21_g2</th><td>  111.314252</td><td>-0.065350469</td><td>0.11493922</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c3_g1</th><td>  211.290828</td><td> 0.192138292</td><td>0.09039701</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c4_g1</th><td>   68.886547</td><td>-0.219215482</td><td>0.21486246</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c4_g2</th><td>22331.591246</td><td> 0.219929202</td><td>0.03636073</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c5_g2</th><td> 2530.149038</td><td> 0.007362662</td><td>0.06621792</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c7_g1</th><td>   79.156674</td><td> 0.241532282</td><td>0.13662389</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c8_g1</th><td> 4418.987496</td><td> 0.134409472</td><td>0.03762453</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN0_c9_g1</th><td>   55.916479</td><td> 0.208580671</td><td>0.12911663</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN1_c1_g2</th><td> 6823.191368</td><td>-0.226669414</td><td>0.07055215</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN1_c2_g1</th><td>   16.643364</td><td>-0.556051373</td><td>0.38818471</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN1_c20_g1</th><td>    5.152759</td><td>-0.259906703</td><td>0.64181022</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN10_c0_g1</th><td> 5907.743404</td><td>-0.010554997</td><td>0.02708673</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c0_g2</th><td> 2100.156725</td><td> 0.128001339</td><td>0.04401182</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c1_g1</th><td>  142.809072</td><td>-0.198911213</td><td>0.08352253</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c1_g3</th><td>    7.738695</td><td> 0.924040144</td><td>0.42304871</td><td> 0.8132400</td><td>0.41608045</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c2_g1</th><td> 1852.416572</td><td>-0.003694662</td><td>0.04707327</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c2_g2</th><td>12806.028272</td><td>-0.045538074</td><td>0.03134368</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c2_g3</th><td>   32.568520</td><td>-0.071390834</td><td>0.20164002</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c4_g1</th><td>    7.025106</td><td> 0.090404364</td><td>0.42233210</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN100_c6_g1</th><td>  350.274326</td><td>-0.471819841</td><td>0.07455599</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN1000_c0_g1</th><td> 3425.272804</td><td> 0.364079705</td><td>0.03313191</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>TRINITY_DN994_c5_g1</th><td>   3.667495</td><td> 0.05007459</td><td>0.45635233</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9941_c0_g1</th><td>  11.475510</td><td> 2.11227073</td><td>0.60992114</td><td> 2.51224400</td><td>0.01199661</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN99413_c0_g1</th><td>  12.710134</td><td> 0.15901563</td><td>0.31886840</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9942_c0_g1</th><td>   7.525042</td><td>-0.62233573</td><td>0.45227485</td><td>-0.09360619</td><td>0.92542199</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN99429_c0_g1</th><td>  90.545316</td><td>-0.04628323</td><td>0.13140714</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9943_c0_g1</th><td>  23.639197</td><td> 0.55942322</td><td>0.20922604</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9946_c0_g2</th><td>  10.267645</td><td> 0.15027193</td><td>0.39868859</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9949_c0_g1</th><td>  12.921730</td><td>-0.36038790</td><td>0.45263255</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN995_c0_g1</th><td> 204.814695</td><td>-0.19168843</td><td>0.09196291</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN995_c1_g1</th><td>   5.118045</td><td> 1.36122006</td><td>0.67251967</td><td> 1.16163153</td><td>0.24538517</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9952_c0_g1</th><td> 285.904262</td><td> 0.12306739</td><td>0.15787162</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9953_c0_g2</th><td>  88.264371</td><td> 0.09916924</td><td>0.18008228</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN996_c1_g1</th><td>  40.733987</td><td> 0.17545432</td><td>0.47467798</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9962_c0_g3</th><td>  13.595194</td><td> 0.48665882</td><td>0.32965248</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN997_c0_g3</th><td>  21.848769</td><td>-0.48491941</td><td>0.28825167</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN99759_c0_g1</th><td>  92.053494</td><td> 0.20530985</td><td>0.11842520</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN99790_c0_g1</th><td>   9.098529</td><td> 0.33553979</td><td>0.35024460</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN998_c0_g1</th><td> 272.405776</td><td>-0.10482563</td><td>0.07448690</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN998_c1_g1</th><td>2507.792218</td><td> 0.07044966</td><td>0.03003691</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9980_c0_g1</th><td>  99.680791</td><td> 0.15200640</td><td>0.09706556</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9984_c0_g1</th><td> 189.712394</td><td> 0.15984120</td><td>0.09915854</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9986_c1_g1</th><td>   4.793399</td><td> 0.24424601</td><td>0.59766435</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN999_c0_g1</th><td>  13.147643</td><td>-0.06295881</td><td>0.28812429</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN99904_c1_g1</th><td>   6.671768</td><td> 0.14528915</td><td>0.60054748</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9991_c0_g1</th><td>  57.855916</td><td> 0.12116419</td><td>0.15767282</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9993_c0_g2</th><td>  77.248675</td><td>-0.11538146</td><td>0.17962431</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9994_c3_g1</th><td>  43.071385</td><td>-0.07261691</td><td>0.19921284</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9995_c0_g1</th><td>   7.219157</td><td> 0.38877925</td><td>0.48482995</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9999_c0_g1</th><td>  37.655719</td><td> 0.11657106</td><td>0.18765409</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><th scope=row>TRINITY_DN9999_c1_g1</th><td>   9.670259</td><td>-0.76154938</td><td>0.38178856</td><td>-0.47552335</td><td>0.63441399</td><td>1</td></tr>
</tbody>
</table>




```R
plotMA(res_table_1, ylim=c(-3,3))
```


    
![png](output_28_0.png)
    



```R
plotMA(res_table_8, ylim=c(-3,3))
```


    
![png](output_29_0.png)
    



```R
plotMA(res_table_24, ylim=c(-3,3))
```


    
![png](output_30_0.png)
    



```R
plotMA(res_table_96, ylim=c(-3,3))
```


    
![png](output_31_0.png)
    



```R
summary(res_table_1, alpha = 0.05)
```

    
    out of 17779 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0.58 (up)    : 53, 0.3%
    LFC < -0.58 (down) : 49, 0.28%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 2)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    


```R
summary(res_table_8, alpha = 0.05)
```

    
    out of 17779 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0.58 (up)    : 0, 0%
    LFC < -0.58 (down) : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 2)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    


```R
summary(res_table_24, alpha = 0.05)
```

    
    out of 17779 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0.58 (up)    : 0, 0%
    LFC < -0.58 (down) : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 2)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    


```R
summary(res_table_96, alpha = 0.05)
```

    
    out of 17779 with nonzero total read count
    adjusted p-value < 0.05
    LFC > 0.58 (up)    : 0, 0%
    LFC < -0.58 (down) : 0, 0%
    outliers [1]       : 0, 0%
    low counts [2]     : 0, 0%
    (mean count < 2)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results
    
    


```R
### Set thresholds
padj.cutoff <- 0.05
```


```R
res_table_1_tb <- res_table_1 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
res_table_8_tb <- res_table_8 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
res_table_24_tb <- res_table_24 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
res_table_96_tb <- res_table_96 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```


```R
sig_1 <- res_table_1_tb %>%
        filter(padj < padj.cutoff)
```


```R
sig_8 <- res_table_8_tb %>%
        filter(padj < padj.cutoff)
```


```R
sig_24 <- res_table_24_tb %>%
        filter(padj < padj.cutoff)
```


```R
sig_96 <- res_table_96_tb %>%
        filter(padj < padj.cutoff)
```


```R
meta_tb <- meta %>% 
              rownames_to_column(var="samplename") %>% 
              as_tibble()
```


```R
DEGreport::degPlot(dds = dds, res = sig_1, genes = rownames(dds)[1:10], xs = "group") # dds object is output from DESeq2
```

    No genes were mapped to rowData. check ann parameter values.
    
    Using gene as id variables
    
    


    Error in `$<-.data.frame`(`*tmp*`, "xs", value = structure(integer(0), .Label = character(0), class = "factor")): replacement has 0 rows, data has 480
    Traceback:
    

    1. DEGreport::degPlot(dds = dds, res = sig_1, genes = rownames(dds)[1:10], 
     .     xs = "group")

    2. `$<-`(`*tmp*`, "xs", value = structure(integer(0), .Label = character(0), class = "factor"))

    3. `$<-.data.frame`(`*tmp*`, "xs", value = structure(integer(0), .Label = character(0), class = "factor"))

    4. stop(sprintf(ngettext(N, "replacement has %d row, data has %d", 
     .     "replacement has %d rows, data has %d"), N, nrows), domain = NA)



```R

```
