```R
#BiocManager::install("Biobase")
#install.packages(c('ggplot2','gplots','NMF','ggfortify'))
#install.packages('ggbreak')
#BiocManager::install("graph")
#library(devtools)
#install_github("igraph/irgraph")
#BiocManager::install("clusterProfiler", force = TRUE)
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
#library(clusterProfiler)

library(ggplot2)
library(gplots)
library(stats)
library(NMF)
library(ggfortify)
library(ggbreak)
library(dplyr)
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
    
    
    Attaching package: 'gplots'
    
    
    The following object is masked from 'package:IRanges':
    
        space
    
    
    The following object is masked from 'package:S4Vectors':
    
        space
    
    
    The following object is masked from 'package:stats':
    
        lowess
    
    
    Loading required package: registry
    
    Loading required package: rngtools
    
    Loading required package: cluster
    
    NMF - BioConductor layer [OK] | Shared memory capabilities [NO: windows] | Cores 2/2
    
    
    Attaching package: 'NMF'
    
    
    The following object is masked from 'package:S4Vectors':
    
        nrun
    
    
    Registered S3 method overwritten by 'ggfortify':
      method        from 
      fortify.table ggalt
    
    [90mggbreak v0.1.2
    
    If you use ggbreak in published research, please cite the following
    paper:
    
    S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively
    utilize plotting space to deal with large datasets and outliers.
    Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846
    [39m
    
    


```R
R.Version()
```


<dl>
	<dt>$platform</dt>
		<dd>'x86_64-w64-mingw32'</dd>
	<dt>$arch</dt>
		<dd>'x86_64'</dd>
	<dt>$os</dt>
		<dd>'mingw32'</dd>
	<dt>$system</dt>
		<dd>'x86_64, mingw32'</dd>
	<dt>$status</dt>
		<dd>''</dd>
	<dt>$major</dt>
		<dd>'4'</dd>
	<dt>$minor</dt>
		<dd>'1.3'</dd>
	<dt>$year</dt>
		<dd>'2022'</dd>
	<dt>$month</dt>
		<dd>'03'</dd>
	<dt>$day</dt>
		<dd>'10'</dd>
	<dt>$`svn rev`</dt>
		<dd>'81868'</dd>
	<dt>$language</dt>
		<dd>'R'</dd>
	<dt>$version.string</dt>
		<dd>'R version 4.1.3 (2022-03-10)'</dd>
	<dt>$nickname</dt>
		<dd>'One Push-Up'</dd>
</dl>




```R
result_table_1 <- read.csv(file = "deseq2_results/results_timepoint_1_swissnames.csv", header = TRUE)
result_table_1
```

    Warning message in file(file, "rt"):
    "cannot open file 'deseq2_results/results_timepoint_1_swissnames.csv': No such file or directory"
    


    Error in file(file, "rt"): cannot open the connection
    Traceback:
    

    1. read.csv(file = "deseq2_results/results_timepoint_1_swissnames.csv", 
     .     header = TRUE)

    2. read.table(file = file, header = header, sep = sep, quote = quote, 
     .     dec = dec, fill = fill, comment.char = comment.char, ...)

    3. file(file, "rt")



```R
# Dropping X Col
result_table_1 <- result_table_1 %>%
  select(-X)

# making trinity gene id col and reording
#result_table_1$trinity_gene <- sapply(strsplit(result_table_1$transcriptid, "_"), function(x) paste(x[1:3], collapse = "_"))

#result_table_1 <- result_table_1 %>%
#  select(trinity_gene, everything())
                            
```


```R
result_table_1
```


<table class="dataframe">
<caption>A data.frame: 17779 × 9</caption>
<thead>
	<tr><th scope=col>transcriptid</th><th scope=col>gene</th><th scope=col>gene_simple</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TRINITY_DN0_c0_g1   </td><td>Unidentified</td><td>Unidentified</td><td> 4898.712104</td><td> 0.210116016</td><td>0.02331269</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c0_g2   </td><td>Unidentified</td><td>Unidentified</td><td>    9.299546</td><td> 0.438524391</td><td>0.34326004</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c1_g1   </td><td>Unidentified</td><td>Unidentified</td><td>  259.542125</td><td>-0.367983662</td><td>0.10516348</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c10_g1  </td><td>RDGB_DROME  </td><td>RDGB        </td><td> 1446.017687</td><td>-0.166869963</td><td>0.03026031</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c11_g1  </td><td>Unidentified</td><td>Unidentified</td><td>    6.576095</td><td> 0.753797609</td><td>0.41938029</td><td> 0.4144153</td><td>0.67856998</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c142_g1 </td><td>PGBM_HUMAN  </td><td>PGBM        </td><td>   84.534552</td><td> 0.196676960</td><td>0.12715418</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c15_g1  </td><td>Unidentified</td><td>Unidentified</td><td>   71.669703</td><td> 0.011562284</td><td>0.14955774</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c2_g1   </td><td>Unidentified</td><td>Unidentified</td><td>    5.545240</td><td>-1.530496463</td><td>0.48734815</td><td>-1.9503439</td><td>0.05113515</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c20_g1  </td><td>Unidentified</td><td>Unidentified</td><td>   35.304042</td><td> 0.537921960</td><td>0.20142297</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c21_g2  </td><td>Unidentified</td><td>Unidentified</td><td>  111.314252</td><td>-0.065350469</td><td>0.11493922</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c3_g1   </td><td>Unidentified</td><td>Unidentified</td><td>  211.290828</td><td> 0.192138292</td><td>0.09039701</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c4_g1   </td><td>Unidentified</td><td>Unidentified</td><td>   68.886547</td><td>-0.219215482</td><td>0.21486246</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c4_g2   </td><td>ABCA3_HUMAN </td><td>ABCA3       </td><td>22331.591246</td><td> 0.219929202</td><td>0.03636073</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c5_g2   </td><td>Unidentified</td><td>Unidentified</td><td> 2530.149038</td><td> 0.007362662</td><td>0.06621792</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c7_g1   </td><td>Unidentified</td><td>Unidentified</td><td>   79.156674</td><td> 0.241532282</td><td>0.13662389</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c8_g1   </td><td>NBEA_DROME  </td><td>NBEA        </td><td> 4418.987496</td><td> 0.134409472</td><td>0.03762453</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN0_c9_g1   </td><td>Unidentified</td><td>Unidentified</td><td>   55.916479</td><td> 0.208580671</td><td>0.12911663</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN1_c1_g2   </td><td>Unidentified</td><td>Unidentified</td><td> 6823.191368</td><td>-0.226669414</td><td>0.07055215</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN1_c2_g1   </td><td>Unidentified</td><td>Unidentified</td><td>   16.643364</td><td>-0.556051373</td><td>0.38818471</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN1_c20_g1  </td><td>Unidentified</td><td>Unidentified</td><td>    5.152759</td><td>-0.259906703</td><td>0.64181022</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN10_c0_g1  </td><td>NPL4_DROME  </td><td>NPL4        </td><td> 5907.743404</td><td>-0.010554997</td><td>0.02708673</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c0_g2 </td><td>RTXE_DROME  </td><td>RTXE        </td><td> 2100.156725</td><td> 0.128001339</td><td>0.04401182</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c1_g1 </td><td>CHN_DROME   </td><td>CHN         </td><td>  142.809072</td><td>-0.198911213</td><td>0.08352253</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c1_g3 </td><td>EGR1_DANRE  </td><td>EGR1        </td><td>    7.738695</td><td> 0.924040144</td><td>0.42304871</td><td> 0.8132400</td><td>0.41608045</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c2_g1 </td><td>Unidentified</td><td>Unidentified</td><td> 1852.416572</td><td>-0.003694662</td><td>0.04707327</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c2_g2 </td><td>TCPH_CHICK  </td><td>TCPH        </td><td>12806.028272</td><td>-0.045538074</td><td>0.03134368</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c2_g3 </td><td>TRYP_ASTAS  </td><td>TRYP        </td><td>   32.568520</td><td>-0.071390834</td><td>0.20164002</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c4_g1 </td><td>Unidentified</td><td>Unidentified</td><td>    7.025106</td><td> 0.090404364</td><td>0.42233210</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN100_c6_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  350.274326</td><td>-0.471819841</td><td>0.07455599</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN1000_c0_g1</td><td>FA13A_MOUSE </td><td>FA13A       </td><td> 3425.272804</td><td> 0.364079705</td><td>0.03313191</td><td> 0.0000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>TRINITY_DN994_c5_g1  </td><td>Unidentified</td><td>Unidentified</td><td>   3.667495</td><td> 0.05007459</td><td>0.45635233</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9941_c0_g1 </td><td>CBPA1_ANOGA </td><td>CBPA1       </td><td>  11.475510</td><td> 2.11227073</td><td>0.60992114</td><td> 2.51224400</td><td>0.01199661</td><td>1</td></tr>
	<tr><td>TRINITY_DN99413_c0_g1</td><td>Unidentified</td><td>Unidentified</td><td>  12.710134</td><td> 0.15901563</td><td>0.31886840</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9942_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td>   7.525042</td><td>-0.62233573</td><td>0.45227485</td><td>-0.09360619</td><td>0.92542199</td><td>1</td></tr>
	<tr><td>TRINITY_DN99429_c0_g1</td><td>Unidentified</td><td>Unidentified</td><td>  90.545316</td><td>-0.04628323</td><td>0.13140714</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9943_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  23.639197</td><td> 0.55942322</td><td>0.20922604</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9946_c0_g2 </td><td>Unidentified</td><td>Unidentified</td><td>  10.267645</td><td> 0.15027193</td><td>0.39868859</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9949_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  12.921730</td><td>-0.36038790</td><td>0.45263255</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN995_c0_g1  </td><td>Unidentified</td><td>Unidentified</td><td> 204.814695</td><td>-0.19168843</td><td>0.09196291</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN995_c1_g1  </td><td>Unidentified</td><td>Unidentified</td><td>   5.118045</td><td> 1.36122006</td><td>0.67251967</td><td> 1.16163153</td><td>0.24538517</td><td>1</td></tr>
	<tr><td>TRINITY_DN9952_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td> 285.904262</td><td> 0.12306739</td><td>0.15787162</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9953_c0_g2 </td><td>Unidentified</td><td>Unidentified</td><td>  88.264371</td><td> 0.09916924</td><td>0.18008228</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN996_c1_g1  </td><td>CUD4_LOCMI  </td><td>CUD4        </td><td>  40.733987</td><td> 0.17545432</td><td>0.47467798</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9962_c0_g3 </td><td>Unidentified</td><td>Unidentified</td><td>  13.595194</td><td> 0.48665882</td><td>0.32965248</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN997_c0_g3  </td><td>Unidentified</td><td>Unidentified</td><td>  21.848769</td><td>-0.48491941</td><td>0.28825167</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN99759_c0_g1</td><td>Unidentified</td><td>Unidentified</td><td>  92.053494</td><td> 0.20530985</td><td>0.11842520</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN99790_c0_g1</td><td>SYTC_HUMAN  </td><td>SYTC        </td><td>   9.098529</td><td> 0.33553979</td><td>0.35024460</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN998_c0_g1  </td><td>Unidentified</td><td>Unidentified</td><td> 272.405776</td><td>-0.10482563</td><td>0.07448690</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN998_c1_g1  </td><td>AXN_DROME   </td><td>AXN         </td><td>2507.792218</td><td> 0.07044966</td><td>0.03003691</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9980_c0_g1 </td><td>PO210_DROME </td><td>PO210       </td><td>  99.680791</td><td> 0.15200640</td><td>0.09706556</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9984_c0_g1 </td><td>CASP6_CHICK </td><td>CASP6       </td><td> 189.712394</td><td> 0.15984120</td><td>0.09915854</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9986_c1_g1 </td><td>Unidentified</td><td>Unidentified</td><td>   4.793399</td><td> 0.24424601</td><td>0.59766435</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN999_c0_g1  </td><td>Unidentified</td><td>Unidentified</td><td>  13.147643</td><td>-0.06295881</td><td>0.28812429</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN99904_c1_g1</td><td>Unidentified</td><td>Unidentified</td><td>   6.671768</td><td> 0.14528915</td><td>0.60054748</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9991_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  57.855916</td><td> 0.12116419</td><td>0.15767282</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9993_c0_g2 </td><td>Unidentified</td><td>Unidentified</td><td>  77.248675</td><td>-0.11538146</td><td>0.17962431</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9994_c3_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  43.071385</td><td>-0.07261691</td><td>0.19921284</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9995_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td>   7.219157</td><td> 0.38877925</td><td>0.48482995</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9999_c0_g1 </td><td>TDC1_CAEEL  </td><td>TDC1        </td><td>  37.655719</td><td> 0.11657106</td><td>0.18765409</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>TRINITY_DN9999_c1_g1 </td><td>Unidentified</td><td>Unidentified</td><td>   9.670259</td><td>-0.76154938</td><td>0.38178856</td><td>-0.47552335</td><td>0.63441399</td><td>1</td></tr>
</tbody>
</table>




```R
norm_count <- read.csv(file = "deseq2_results/norm_counts/Normalized_Counts_avgTxLength.csv", header = TRUE)
norm_count <- norm_count %>%
  rename(transcriptid = X)
```


```R
norm_count
```


<table class="dataframe">
<caption>A data.frame: 17779 × 49</caption>
<thead>
	<tr><th scope=col>transcriptid</th><th scope=col>C_1_A</th><th scope=col>C_1_B</th><th scope=col>C_1_C</th><th scope=col>C_1_D</th><th scope=col>C_1_E</th><th scope=col>C_1_F</th><th scope=col>C_8_A</th><th scope=col>C_8_B</th><th scope=col>C_8_C</th><th scope=col>...</th><th scope=col>V_24_C</th><th scope=col>V_24_D</th><th scope=col>V_24_E</th><th scope=col>V_24_F</th><th scope=col>V_96_A</th><th scope=col>V_96_B</th><th scope=col>V_96_C</th><th scope=col>V_96_D</th><th scope=col>V_96_E</th><th scope=col>V_96_F</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>...</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TRINITY_DN0_c0_g1   </td><td>2996.41886</td><td>2985.34643</td><td>3023.59815</td><td>2921.25842</td><td>2945.80908</td><td>3000.1189</td><td>2745.60207</td><td>2768.28017</td><td>2810.53457</td><td>...</td><td>2996.16462</td><td>2989.97606</td><td>2981.66212</td><td>3021.60509</td><td>2662.42684</td><td>2638.14957</td><td>2707.64967</td><td>2718.08463</td><td>2727.83861</td><td>2696.93756</td></tr>
	<tr><td>TRINITY_DN0_c0_g2   </td><td> 284.09500</td><td> 298.82000</td><td> 307.71600</td><td> 287.63300</td><td> 295.19400</td><td> 293.6300</td><td> 301.90100</td><td> 301.02100</td><td> 314.36700</td><td>...</td><td> 294.22900</td><td> 292.78400</td><td> 295.87300</td><td> 290.44200</td><td> 299.24100</td><td> 297.42700</td><td> 299.79000</td><td> 305.38200</td><td> 299.91500</td><td> 294.15200</td></tr>
	<tr><td>TRINITY_DN0_c1_g1   </td><td> 649.33900</td><td> 664.27800</td><td> 673.27700</td><td> 653.13900</td><td> 660.70100</td><td> 659.1350</td><td> 667.38500</td><td> 666.48200</td><td> 679.94800</td><td>...</td><td> 659.60900</td><td> 658.23600</td><td> 661.35200</td><td> 655.87700</td><td> 664.72700</td><td> 662.83300</td><td> 665.13800</td><td> 670.88500</td><td> 665.37500</td><td> 659.57600</td></tr>
	<tr><td>TRINITY_DN0_c10_g1  </td><td>2523.31800</td><td>2538.26100</td><td>2547.25800</td><td>2527.12300</td><td>2534.68100</td><td>2533.1160</td><td>2541.35600</td><td>2540.46400</td><td>2553.93200</td><td>...</td><td>2533.58900</td><td>2532.21900</td><td>2535.33100</td><td>2529.85600</td><td>2538.70400</td><td>2536.81200</td><td>2539.12500</td><td>2544.86100</td><td>2539.35500</td><td>2533.56300</td></tr>
	<tr><td>TRINITY_DN0_c11_g1  </td><td>  99.78500</td><td> 109.77600</td><td> 116.39900</td><td> 102.71800</td><td> 106.71100</td><td> 106.1220</td><td> 113.70200</td><td> 111.88200</td><td> 121.81200</td><td>...</td><td> 106.05300</td><td> 105.38800</td><td> 107.28800</td><td> 103.47200</td><td> 107.55600</td><td> 106.34000</td><td> 108.54500</td><td> 112.89400</td><td> 108.20600</td><td> 104.57700</td></tr>
	<tr><td>TRINITY_DN0_c142_g1 </td><td> 226.67500</td><td> 241.22700</td><td> 250.01200</td><td> 230.13200</td><td> 237.57000</td><td> 236.0430</td><td> 244.28600</td><td> 243.41300</td><td> 256.62500</td><td>...</td><td> 236.66700</td><td> 235.24000</td><td> 238.28100</td><td> 232.90800</td><td> 241.55600</td><td> 239.80800</td><td> 242.18800</td><td> 247.68300</td><td> 242.23000</td><td> 236.55800</td></tr>
	<tr><td>TRINITY_DN0_c15_g1  </td><td> 204.15500</td><td> 218.51900</td><td> 227.21600</td><td> 207.54500</td><td> 214.85500</td><td> 213.3710</td><td> 221.61000</td><td> 220.72500</td><td> 233.80200</td><td>...</td><td> 214.01700</td><td> 212.58800</td><td> 215.59700</td><td> 210.29500</td><td> 218.78300</td><td> 217.06100</td><td> 219.45200</td><td> 224.88900</td><td> 219.46600</td><td> 213.85000</td></tr>
	<tr><td>TRINITY_DN0_c2_g1   </td><td>  58.53600</td><td>  59.52900</td><td>  61.68000</td><td>  61.26900</td><td>  58.08100</td><td>  60.6630</td><td>  66.73500</td><td>  61.76000</td><td>  64.70500</td><td>...</td><td>  57.44000</td><td>  59.13700</td><td>  58.74200</td><td>  58.18100</td><td>  51.64200</td><td>  50.91100</td><td>  58.30643</td><td>  54.42600</td><td>  52.18400</td><td>  53.14600</td></tr>
	<tr><td>TRINITY_DN0_c20_g1  </td><td>  59.58400</td><td>  61.00900</td><td>  63.34600</td><td>  62.36500</td><td>  59.46500</td><td>  61.9450</td><td>  68.02400</td><td>  63.18500</td><td>  66.45800</td><td>...</td><td>  58.80200</td><td>  60.37600</td><td>  60.05700</td><td>  59.35700</td><td>  53.36600</td><td>  52.62300</td><td>  54.07700</td><td>  56.29600</td><td>  53.95200</td><td>  54.63800</td></tr>
	<tr><td>TRINITY_DN0_c21_g2  </td><td> 136.63500</td><td> 149.37500</td><td> 157.37600</td><td> 139.77500</td><td> 145.93000</td><td> 144.7090</td><td> 152.71800</td><td> 151.57300</td><td> 163.57100</td><td>...</td><td> 145.25300</td><td> 144.00700</td><td> 146.65900</td><td> 141.86400</td><td> 148.85500</td><td> 147.32700</td><td> 149.71700</td><td> 154.73300</td><td> 149.57000</td><td> 144.57600</td></tr>
	<tr><td>TRINITY_DN0_c3_g1   </td><td>  89.31773</td><td>  91.18942</td><td>  99.16362</td><td>  95.17226</td><td>  96.08242</td><td>  92.4825</td><td>  96.55127</td><td>  94.14659</td><td>  90.39355</td><td>...</td><td>  86.49768</td><td>  84.07179</td><td>  89.28264</td><td>  81.55056</td><td>  78.00351</td><td>  74.98678</td><td>  78.55074</td><td>  83.33545</td><td>  82.21145</td><td>  80.32738</td></tr>
	<tr><td>TRINITY_DN0_c4_g1   </td><td>3817.60722</td><td>4039.45190</td><td>4336.55590</td><td>3753.28776</td><td>5404.76055</td><td>3696.6130</td><td>3311.55181</td><td>3456.74306</td><td>3006.48122</td><td>...</td><td>3369.05774</td><td>3017.10274</td><td>3435.75784</td><td>3292.30074</td><td>6534.26532</td><td>5865.75242</td><td>5629.61432</td><td>6659.42781</td><td>6841.31457</td><td>6538.12806</td></tr>
	<tr><td>TRINITY_DN0_c4_g2   </td><td>5566.31800</td><td>5581.26100</td><td>5590.25800</td><td>5570.12300</td><td>5577.68100</td><td>5576.1160</td><td>5584.35600</td><td>5583.46400</td><td>5596.93200</td><td>...</td><td>5576.58900</td><td>5575.21900</td><td>5578.33100</td><td>5572.85600</td><td>5581.70400</td><td>5579.81200</td><td>5582.12500</td><td>5587.86100</td><td>5582.35500</td><td>5576.56300</td></tr>
	<tr><td>TRINITY_DN0_c5_g2   </td><td> 121.31400</td><td> 133.17300</td><td> 140.77300</td><td> 124.36000</td><td> 129.84100</td><td> 128.8020</td><td> 136.65200</td><td> 135.34600</td><td> 146.72600</td><td>...</td><td> 129.20000</td><td> 128.13100</td><td> 130.54900</td><td> 126.04700</td><td> 132.14600</td><td> 130.73200</td><td> 133.05700</td><td> 137.86700</td><td> 132.84100</td><td> 128.25600</td></tr>
	<tr><td>TRINITY_DN0_c7_g1   </td><td> 320.88700</td><td> 335.69400</td><td> 344.62400</td><td> 324.49000</td><td> 332.06500</td><td> 330.4990</td><td> 338.77100</td><td> 337.89400</td><td> 351.27000</td><td>...</td><td> 331.07700</td><td> 329.63900</td><td> 332.75000</td><td> 327.29700</td><td> 336.11200</td><td> 334.28400</td><td> 336.62200</td><td> 342.26300</td><td> 336.78600</td><td> 331.00100</td></tr>
	<tr><td>TRINITY_DN0_c8_g1   </td><td>3917.31800</td><td>3932.26100</td><td>3941.25800</td><td>3921.12300</td><td>3928.68100</td><td>3927.1160</td><td>3935.35600</td><td>3934.46400</td><td>3947.93200</td><td>...</td><td>3927.58900</td><td>3926.21900</td><td>3929.33100</td><td>3923.85600</td><td>3932.70400</td><td>3930.81200</td><td>3933.12500</td><td>3938.86100</td><td>3933.35500</td><td>3927.56300</td></tr>
	<tr><td>TRINITY_DN0_c9_g1   </td><td> 477.54300</td><td> 492.42300</td><td> 501.39800</td><td> 481.28400</td><td> 488.83400</td><td> 487.2780</td><td> 495.53800</td><td> 494.63900</td><td> 508.06300</td><td>...</td><td> 487.77200</td><td> 486.36800</td><td> 489.49400</td><td> 484.04400</td><td> 492.85200</td><td> 490.97100</td><td> 493.28500</td><td> 499.00600</td><td> 493.52100</td><td> 487.71900</td></tr>
	<tr><td>TRINITY_DN1_c1_g2   </td><td>1877.50736</td><td>1892.08156</td><td>1900.30884</td><td>1881.20741</td><td>1888.45447</td><td>1887.5068</td><td>1897.66166</td><td>1896.24409</td><td>1908.48293</td><td>...</td><td>1888.79609</td><td>1887.02838</td><td>1890.30027</td><td>1885.01644</td><td>1888.91155</td><td>1886.91741</td><td>1889.40348</td><td>1895.09519</td><td>1889.53236</td><td>1883.68383</td></tr>
	<tr><td>TRINITY_DN1_c2_g1   </td><td> 128.98629</td><td>  89.23869</td><td> 151.89092</td><td>  81.72907</td><td>  97.67150</td><td> 146.6085</td><td> 110.45585</td><td>  96.70345</td><td> 102.24748</td><td>...</td><td> 150.41165</td><td> 154.73262</td><td> 136.29158</td><td> 100.32467</td><td> 331.46381</td><td> 175.97990</td><td> 197.23720</td><td> 751.86100</td><td> 171.47717</td><td> 167.20697</td></tr>
	<tr><td>TRINITY_DN1_c20_g1  </td><td> 974.41552</td><td> 975.26100</td><td> 984.25800</td><td> 964.12300</td><td> 971.68100</td><td> 970.1160</td><td> 974.41552</td><td> 974.41552</td><td> 990.93200</td><td>...</td><td> 970.58900</td><td> 969.21900</td><td> 974.41552</td><td> 966.85600</td><td> 975.70400</td><td> 973.81200</td><td> 976.12500</td><td> 981.86100</td><td> 976.35500</td><td> 970.56300</td></tr>
	<tr><td>TRINITY_DN10_c0_g1  </td><td>2995.72149</td><td>3025.00585</td><td>3025.98928</td><td>3010.95056</td><td>3023.08895</td><td>3005.0750</td><td>3126.91148</td><td>3153.06474</td><td>3133.84127</td><td>...</td><td>2992.64479</td><td>2993.33970</td><td>3014.85807</td><td>2996.28513</td><td>3183.17670</td><td>3171.98911</td><td>3234.87021</td><td>3215.08753</td><td>3161.81146</td><td>3193.54949</td></tr>
	<tr><td>TRINITY_DN100_c0_g2 </td><td>4174.61509</td><td>4186.72920</td><td>4202.75313</td><td>4174.62089</td><td>4190.47582</td><td>4188.5693</td><td>4184.34641</td><td>4200.45153</td><td>4210.86130</td><td>...</td><td>4180.51766</td><td>4181.45411</td><td>4180.37953</td><td>4173.80371</td><td>4195.46853</td><td>4191.03710</td><td>4191.17902</td><td>4197.60636</td><td>4196.73767</td><td>4186.50363</td></tr>
	<tr><td>TRINITY_DN100_c1_g1 </td><td>1287.31800</td><td>1302.26100</td><td>1311.25800</td><td>1291.12300</td><td>1298.68100</td><td>1297.1160</td><td>1305.35600</td><td>1304.46400</td><td>1317.93200</td><td>...</td><td>1297.58900</td><td>1296.21900</td><td>1299.33100</td><td>1293.85600</td><td>1302.70400</td><td>1300.81200</td><td>1303.12500</td><td>1308.86100</td><td>1303.35500</td><td>1297.56300</td></tr>
	<tr><td>TRINITY_DN100_c1_g3 </td><td>2258.31800</td><td>2273.26100</td><td>2282.25800</td><td>2262.12300</td><td>2269.68100</td><td>2268.1160</td><td>2276.35600</td><td>2275.46400</td><td>2288.93200</td><td>...</td><td>2268.58900</td><td>2267.21900</td><td>2270.33100</td><td>2264.85600</td><td>2273.70400</td><td>2271.81200</td><td>2274.12500</td><td>2279.86100</td><td>2274.35500</td><td>2268.56300</td></tr>
	<tr><td>TRINITY_DN100_c2_g1 </td><td>1939.19853</td><td>1997.84213</td><td>1877.12528</td><td>1762.52902</td><td>1811.91096</td><td>1848.3369</td><td>2186.58328</td><td>1962.05074</td><td>2171.26344</td><td>...</td><td>2221.48027</td><td>2190.71332</td><td>2081.45450</td><td>1967.45275</td><td>1976.24061</td><td>1997.96088</td><td>1972.17847</td><td>2008.32483</td><td>1965.07894</td><td>1939.12285</td></tr>
	<tr><td>TRINITY_DN100_c2_g2 </td><td>2266.31800</td><td>2281.26100</td><td>2290.25800</td><td>2270.12300</td><td>2277.68100</td><td>2276.1160</td><td>2284.35600</td><td>2283.46400</td><td>2296.93200</td><td>...</td><td>2276.58900</td><td>2275.21900</td><td>2278.33100</td><td>2272.85600</td><td>2281.70400</td><td>2279.81200</td><td>2282.12500</td><td>2287.86100</td><td>2282.35500</td><td>2276.56300</td></tr>
	<tr><td>TRINITY_DN100_c2_g3 </td><td>1433.31800</td><td>1448.26100</td><td>1457.25800</td><td>1437.12300</td><td>1444.68100</td><td>1443.1160</td><td>1451.35600</td><td>1450.46400</td><td>1463.93200</td><td>...</td><td>1443.58900</td><td>1442.21900</td><td>1445.33100</td><td>1439.85600</td><td>1448.70400</td><td>1446.81200</td><td>1449.12500</td><td>1454.86100</td><td>1449.35500</td><td>1443.56300</td></tr>
	<tr><td>TRINITY_DN100_c4_g1 </td><td>  61.20300</td><td>  63.32900</td><td>  65.94900</td><td>  64.02400</td><td>  61.63300</td><td>  63.8820</td><td>  70.02200</td><td>  65.39500</td><td>  69.17900</td><td>...</td><td>  60.96000</td><td>  62.30600</td><td>  62.17000</td><td>  61.19500</td><td>  56.03400</td><td>  55.30900</td><td>  56.76800</td><td>  59.17600</td><td>  56.66100</td><td>  56.99200</td></tr>
	<tr><td>TRINITY_DN100_c6_g1 </td><td>1524.31800</td><td>1539.26100</td><td>1548.25800</td><td>1528.12300</td><td>1535.68100</td><td>1534.1160</td><td>1542.35600</td><td>1541.46400</td><td>1554.93200</td><td>...</td><td>1534.58900</td><td>1533.21900</td><td>1536.33100</td><td>1530.85600</td><td>1539.70400</td><td>1537.81200</td><td>1540.12500</td><td>1545.86100</td><td>1540.35500</td><td>1534.56300</td></tr>
	<tr><td>TRINITY_DN1000_c0_g1</td><td>2823.31800</td><td>2838.26100</td><td>2847.25800</td><td>2827.12300</td><td>2834.68100</td><td>2833.1160</td><td>2841.35600</td><td>2840.46400</td><td>2853.93200</td><td>...</td><td>2833.58900</td><td>2832.21900</td><td>2835.33100</td><td>2829.85600</td><td>2838.70400</td><td>2836.81200</td><td>2839.12500</td><td>2844.86100</td><td>2839.35500</td><td>2833.56300</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td></td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>TRINITY_DN994_c5_g1  </td><td>1058.3180</td><td>1073.2610</td><td>1082.2580</td><td>1062.123</td><td>1069.681</td><td>1068.11600</td><td>1076.3560</td><td>1075.4640</td><td>1088.9320</td><td>...</td><td>1068.5890</td><td>1067.2190</td><td>1072.0249</td><td>1072.0249</td><td>1073.7040</td><td>1071.81200</td><td>1074.1250</td><td>1079.861</td><td>1074.3550</td><td>1068.56300</td></tr>
	<tr><td>TRINITY_DN9941_c0_g1 </td><td>3072.3180</td><td>3087.2610</td><td>3096.2580</td><td>1000.707</td><td>3083.681</td><td> 757.39866</td><td>3090.3560</td><td>3089.4640</td><td>3102.9320</td><td>...</td><td>3082.5890</td><td>3081.2190</td><td> 411.6062</td><td>3078.8560</td><td>1187.5224</td><td>3085.81200</td><td> 587.7670</td><td>3093.861</td><td>3088.3550</td><td>3082.56300</td></tr>
	<tr><td>TRINITY_DN99413_c0_g1</td><td>  55.8800</td><td>  55.7870</td><td>  57.3180</td><td>  58.454</td><td>  54.456</td><td>  57.34100</td><td>  63.4150</td><td>  58.0190</td><td>  60.2140</td><td>...</td><td>  54.0250</td><td>  55.8830</td><td>  55.3220</td><td>  55.2120</td><td>  47.1760</td><td>  46.56600</td><td>  47.8660</td><td>  49.574</td><td>  47.6700</td><td>  49.35600</td></tr>
	<tr><td>TRINITY_DN9942_c0_g1 </td><td> 787.3180</td><td>1526.1746</td><td>1490.9479</td><td>1502.128</td><td> 798.681</td><td>1593.42348</td><td>1298.3560</td><td>1663.4640</td><td>1023.3094</td><td>...</td><td>1656.5890</td><td>1655.2190</td><td>1368.4394</td><td>1652.8560</td><td>1479.6750</td><td>1024.79370</td><td>1447.7334</td><td>1667.861</td><td>1472.4781</td><td>1173.68918</td></tr>
	<tr><td>TRINITY_DN99429_c0_g1</td><td>  54.4660</td><td>  53.9030</td><td>  55.2400</td><td>  56.831</td><td>  52.565</td><td>  55.62400</td><td>  61.6650</td><td>  56.1380</td><td>  58.0300</td><td>...</td><td>  52.3230</td><td>  54.1210</td><td>  53.6440</td><td>  53.6110</td><td>  44.9840</td><td>  44.44300</td><td>  45.7560</td><td>  47.240</td><td>  45.4450</td><td>  47.47600</td></tr>
	<tr><td>TRINITY_DN9943_c0_g1 </td><td> 101.2190</td><td> 111.3590</td><td> 118.0690</td><td> 104.143</td><td> 108.266</td><td> 107.62500</td><td> 115.2430</td><td> 113.4770</td><td> 123.5180</td><td>...</td><td> 107.6050</td><td> 106.9190</td><td> 108.8570</td><td> 104.9910</td><td> 109.2370</td><td> 108.00100</td><td> 110.2260</td><td> 114.615</td><td> 109.9010</td><td> 106.18400</td></tr>
	<tr><td>TRINITY_DN9946_c0_g2 </td><td>1775.3180</td><td>1790.2610</td><td>1799.2580</td><td>1779.123</td><td>1786.681</td><td>1785.11600</td><td>1793.3560</td><td>1792.4640</td><td>1805.9320</td><td>...</td><td>1785.5890</td><td>1784.2190</td><td>1787.3310</td><td>1781.8560</td><td>1790.7040</td><td>1788.81200</td><td>1791.1250</td><td>1796.861</td><td>1791.3550</td><td>1785.56300</td></tr>
	<tr><td>TRINITY_DN9949_c0_g1 </td><td>1067.3180</td><td>1082.2610</td><td>1091.2580</td><td>1071.123</td><td>1078.681</td><td>1077.11600</td><td>1085.3560</td><td>1084.4640</td><td>1097.9320</td><td>...</td><td>1077.5890</td><td>1076.2190</td><td>1079.3310</td><td>1073.8560</td><td>1082.7040</td><td>1080.81200</td><td>1083.1250</td><td>1088.861</td><td>1083.3550</td><td>1077.56300</td></tr>
	<tr><td>TRINITY_DN995_c0_g1  </td><td> 678.3270</td><td> 693.2740</td><td> 702.2730</td><td> 682.133</td><td> 689.691</td><td> 688.12700</td><td> 696.3760</td><td> 695.4750</td><td> 708.9440</td><td>...</td><td> 688.5990</td><td> 687.2260</td><td> 690.3380</td><td> 684.8680</td><td> 693.7200</td><td> 691.82500</td><td> 694.1300</td><td> 699.877</td><td> 694.3640</td><td> 688.57100</td></tr>
	<tr><td>TRINITY_DN995_c1_g1  </td><td>3575.2931</td><td>3575.2931</td><td>3585.2580</td><td>3565.123</td><td>3572.681</td><td>3571.11600</td><td>3575.2931</td><td>3578.4640</td><td>3591.9320</td><td>...</td><td>3571.5890</td><td>3570.2190</td><td>3573.3310</td><td>3567.8560</td><td>3576.7040</td><td>3574.81200</td><td>3577.1250</td><td>3582.861</td><td>3577.3550</td><td>3571.56300</td></tr>
	<tr><td>TRINITY_DN9952_c0_g1 </td><td>1112.3180</td><td>1131.9957</td><td>1137.1757</td><td>1122.155</td><td>1125.570</td><td>1125.42575</td><td>1133.6664</td><td>1131.8265</td><td>1144.7195</td><td>...</td><td>1128.1162</td><td>1129.5300</td><td>1124.3310</td><td>1122.9958</td><td>1129.1498</td><td>1128.46098</td><td>1129.9316</td><td>1135.228</td><td>1130.4443</td><td>1124.45638</td></tr>
	<tr><td>TRINITY_DN9953_c0_g2 </td><td>1391.3180</td><td>1406.2610</td><td>1415.2580</td><td>1395.123</td><td>1402.681</td><td>1401.11600</td><td>1409.3560</td><td>1408.4640</td><td>1421.9320</td><td>...</td><td>1401.5890</td><td>1400.2190</td><td>1403.3310</td><td>1397.8560</td><td>1406.7040</td><td>1404.81200</td><td>1407.1250</td><td>1412.861</td><td>1407.3550</td><td>1401.56300</td></tr>
	<tr><td>TRINITY_DN996_c1_g1  </td><td>1096.3180</td><td>1111.2610</td><td>1120.2580</td><td>1100.123</td><td>1107.681</td><td>1106.11600</td><td>1114.3560</td><td>1113.4640</td><td>1126.9320</td><td>...</td><td>1106.5890</td><td>1105.2190</td><td>1108.3310</td><td>1102.8560</td><td>1111.7040</td><td>1109.81200</td><td>1112.1250</td><td>1117.861</td><td>1112.3550</td><td>1106.56300</td></tr>
	<tr><td>TRINITY_DN9962_c0_g3 </td><td> 101.2190</td><td> 111.3590</td><td> 118.0690</td><td> 104.143</td><td> 108.266</td><td> 107.62500</td><td> 115.2430</td><td> 113.4770</td><td> 123.5180</td><td>...</td><td> 107.6050</td><td> 106.9190</td><td> 108.8570</td><td> 104.9910</td><td> 109.2370</td><td> 108.00100</td><td> 110.2260</td><td> 114.615</td><td> 109.9010</td><td> 106.18400</td></tr>
	<tr><td>TRINITY_DN997_c0_g3  </td><td>4563.3180</td><td>4578.2610</td><td>4587.2580</td><td>4567.123</td><td>4574.681</td><td>4573.11600</td><td>4581.3560</td><td>4580.4640</td><td>4593.9320</td><td>...</td><td>4573.5890</td><td>4572.2190</td><td>4575.3310</td><td>4569.8560</td><td>4578.7040</td><td>4576.81200</td><td>4579.1250</td><td>4584.861</td><td>4579.3550</td><td>4573.56300</td></tr>
	<tr><td>TRINITY_DN99759_c0_g1</td><td>  63.0170</td><td>  65.8840</td><td>  68.8050</td><td>  65.875</td><td>  64.039</td><td>  65.97800</td><td>  72.2710</td><td>  67.9080</td><td>  72.1880</td><td>...</td><td>  63.3580</td><td>  64.5400</td><td>  64.5530</td><td>  63.3560</td><td>  59.0590</td><td>  58.26700</td><td>  59.8040</td><td>  62.364</td><td>  59.6430</td><td>  59.63200</td></tr>
	<tr><td>TRINITY_DN99790_c0_g1</td><td>  78.7700</td><td>  85.6370</td><td>  90.7040</td><td>  81.628</td><td>  83.161</td><td>  83.53000</td><td>  90.4410</td><td>  87.6990</td><td>  95.1730</td><td>...</td><td>  82.3530</td><td>  82.4340</td><td>  83.6000</td><td>  80.8500</td><td>  81.3790</td><td>  84.29772</td><td>  82.2740</td><td>  85.926</td><td>  82.0210</td><td>  84.29772</td></tr>
	<tr><td>TRINITY_DN998_c0_g1  </td><td>3234.3180</td><td>3249.2610</td><td>3258.2580</td><td>3238.123</td><td>3245.681</td><td>3244.11600</td><td>3252.3560</td><td>3251.4640</td><td>3264.9320</td><td>...</td><td>3244.5890</td><td>3243.2190</td><td>3246.3310</td><td>3240.8560</td><td>3249.7040</td><td>3247.81200</td><td>3250.1250</td><td>3255.861</td><td>3250.3550</td><td>3244.56300</td></tr>
	<tr><td>TRINITY_DN998_c1_g1  </td><td>2708.3180</td><td>2723.2610</td><td>2732.2580</td><td>2712.123</td><td>2719.681</td><td>2718.11600</td><td>2726.3560</td><td>2725.4640</td><td>2738.9320</td><td>...</td><td>2718.5890</td><td>2717.2190</td><td>2720.3310</td><td>2714.8560</td><td>2723.7040</td><td>2721.81200</td><td>2724.1250</td><td>2729.861</td><td>2724.3550</td><td>2718.56300</td></tr>
	<tr><td>TRINITY_DN9980_c0_g1 </td><td> 146.3960</td><td> 159.5990</td><td> 167.7490</td><td> 149.612</td><td> 156.068</td><td> 154.76200</td><td> 162.8720</td><td> 161.7830</td><td> 174.0690</td><td>...</td><td> 155.3830</td><td> 154.0850</td><td> 156.8190</td><td> 151.9080</td><td> 159.3030</td><td> 157.75000</td><td> 160.1060</td><td> 165.220</td><td> 159.9960</td><td> 154.82500</td></tr>
	<tr><td>TRINITY_DN9984_c0_g1 </td><td> 194.4400</td><td> 208.6790</td><td> 217.3430</td><td> 197.806</td><td> 205.049</td><td> 203.57200</td><td> 211.8060</td><td> 210.9010</td><td> 223.9110</td><td>...</td><td> 204.2290</td><td> 202.8030</td><td> 205.7870</td><td> 200.5350</td><td> 208.9210</td><td> 207.22400</td><td> 209.6030</td><td> 215.001</td><td> 209.6070</td><td> 204.01900</td></tr>
	<tr><td>TRINITY_DN9986_c1_g1 </td><td> 230.3128</td><td> 231.3500</td><td> 240.0890</td><td> 220.282</td><td> 227.672</td><td> 226.17100</td><td> 234.4050</td><td> 233.5280</td><td> 246.6900</td><td>...</td><td> 226.8080</td><td> 225.3700</td><td> 228.3990</td><td> 223.0500</td><td> 231.6460</td><td> 229.89800</td><td> 230.3128</td><td> 237.761</td><td> 232.3220</td><td> 226.67100</td></tr>
	<tr><td>TRINITY_DN999_c0_g1  </td><td>1251.4310</td><td>1219.3218</td><td>1334.4078</td><td>1239.464</td><td>1301.952</td><td>1303.76533</td><td>1281.5242</td><td>1253.3663</td><td>1391.4357</td><td>...</td><td>1510.8626</td><td>1451.7658</td><td>1494.2750</td><td>1355.0879</td><td>1357.7094</td><td>1263.14661</td><td>1440.7777</td><td>1363.806</td><td>1573.6360</td><td>1341.25760</td></tr>
	<tr><td>TRINITY_DN99904_c1_g1</td><td>  75.2260</td><td>  81.4260</td><td>  86.0860</td><td>  78.083</td><td>  79.021</td><td>  80.17151</td><td>  86.4410</td><td>  83.4300</td><td>  90.3880</td><td>...</td><td>  78.2370</td><td>  78.5120</td><td>  79.5180</td><td>  76.9900</td><td>  76.6900</td><td>  75.81400</td><td>  77.5510</td><td>  81.012</td><td>  77.3340</td><td>  75.55500</td></tr>
	<tr><td>TRINITY_DN9991_c0_g1 </td><td>  83.6700</td><td>  91.4930</td><td>  96.9440</td><td>  86.569</td><td>  88.822</td><td>  88.87400</td><td>  96.0330</td><td>  93.5620</td><td> 101.7030</td><td>...</td><td>  88.0590</td><td>  87.9130</td><td>  89.2960</td><td>  86.2080</td><td>  87.7790</td><td>  86.80100</td><td>  88.7440</td><td>  92.594</td><td>  88.4750</td><td>  85.86100</td></tr>
	<tr><td>TRINITY_DN9993_c0_g2 </td><td> 372.4719</td><td> 424.1194</td><td> 437.6288</td><td> 419.344</td><td> 366.401</td><td> 390.99998</td><td> 400.7033</td><td> 401.5267</td><td> 444.2533</td><td>...</td><td> 314.7681</td><td> 419.8207</td><td> 328.1724</td><td> 287.9926</td><td> 216.5768</td><td> 196.19307</td><td> 211.6808</td><td> 188.734</td><td> 213.7933</td><td> 288.37829</td></tr>
	<tr><td>TRINITY_DN9994_c3_g1 </td><td>  58.1460</td><td>  59.0010</td><td>  61.0460</td><td>  60.885</td><td>  57.578</td><td>  60.18600</td><td>  66.2530</td><td>  61.2210</td><td>  64.0580</td><td>...</td><td>  56.9330</td><td>  58.6690</td><td>  58.2690</td><td>  57.7560</td><td>  50.9860</td><td>  50.27900</td><td>  51.6870</td><td>  53.719</td><td>  51.5110</td><td>  52.59600</td></tr>
	<tr><td>TRINITY_DN9995_c0_g1 </td><td> 947.3180</td><td> 962.2610</td><td> 971.2580</td><td> 951.123</td><td> 958.681</td><td> 957.11600</td><td> 965.3560</td><td> 964.4640</td><td> 977.9320</td><td>...</td><td> 957.5890</td><td> 956.2190</td><td> 959.3310</td><td> 953.8560</td><td> 962.7040</td><td> 960.81200</td><td> 963.1250</td><td> 968.861</td><td> 963.3550</td><td> 957.56300</td></tr>
	<tr><td>TRINITY_DN9999_c0_g1 </td><td>2231.3960</td><td>2273.2268</td><td>2258.3217</td><td>2264.493</td><td>2291.221</td><td>2246.43144</td><td>2181.2114</td><td>2236.1297</td><td>2249.5778</td><td>...</td><td>2114.5601</td><td>2100.8532</td><td>2048.5254</td><td>2075.1476</td><td>2093.9175</td><td>2103.77103</td><td>2130.9588</td><td>2159.336</td><td>2149.6385</td><td>2134.97050</td></tr>
	<tr><td>TRINITY_DN9999_c1_g1 </td><td> 102.6620</td><td> 112.9520</td><td> 119.7520</td><td> 105.592</td><td> 109.840</td><td> 109.15400</td><td> 116.8030</td><td> 115.0870</td><td> 125.2310</td><td>...</td><td> 109.1780</td><td> 108.4680</td><td> 110.4450</td><td> 106.5180</td><td> 110.9360</td><td> 109.68600</td><td> 111.9160</td><td> 116.346</td><td> 111.6080</td><td> 107.80800</td></tr>
</tbody>
</table>




```R
genes_interest <- list('SMAD6', 'LIPT', 'TOLL7', 'SAM11', 'KI26L', 'KEN1',
                       'Traf4', 'ABCA3', 'CHIT1', 'FAT', 'GGT1', 'PPAF3',
                       'REG5', 'DSX', 'CUD2', 'SCX', 'CP4CU', 'NAAT1',
                       'TRAF4', 'ATG2B', 'NOC', 'GROU', 'USH', 'RIM2',
                       'SCPDL')
```


```R
## Gene of intrest
intrest_sig <- result_table_1 %>%
        filter(gene_simple %in% genes_interest)
```


```R
intrest_sig
```


<table class="dataframe">
<caption>A data.frame: 61 × 9</caption>
<thead>
	<tr><th scope=col>transcriptid</th><th scope=col>gene</th><th scope=col>gene_simple</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TRINITY_DN0_c4_g2    </td><td>ABCA3_HUMAN</td><td>ABCA3</td><td>22331.591246</td><td> 0.21992920</td><td>0.03636073</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN10814_c4_g1</td><td>LIPT_BOVIN </td><td>LIPT </td><td>  286.118561</td><td>-1.34701357</td><td>0.08709560</td><td> -8.806571</td><td>1.290326e-18</td><td>1.147035e-15</td></tr>
	<tr><td>TRINITY_DN1083_c3_g1 </td><td>NAAT1_DROER</td><td>NAAT1</td><td>  103.837277</td><td>-0.25361188</td><td>0.11420024</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN11417_c0_g1</td><td>FAT_DROME  </td><td>FAT  </td><td>   29.765931</td><td> 1.78457691</td><td>0.21336843</td><td>  5.645526</td><td>1.646774e-08</td><td>6.506220e-06</td></tr>
	<tr><td>TRINITY_DN1151_c0_g1 </td><td>PPAF3_HOLDI</td><td>PPAF3</td><td>  123.734742</td><td> 1.29990460</td><td>0.11487874</td><td>  6.266648</td><td>3.689039e-10</td><td>1.725985e-07</td></tr>
	<tr><td>TRINITY_DN1178_c0_g1 </td><td>CHIT1_HUMAN</td><td>CHIT1</td><td> 9207.152270</td><td>-0.42350501</td><td>0.03463782</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN119_c79_g1 </td><td>ABCA3_MOUSE</td><td>ABCA3</td><td> 1946.615355</td><td> 0.26539607</td><td>0.04486380</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN1244_c3_g1 </td><td>LIPT_MOUSE </td><td>LIPT </td><td>21101.545911</td><td>-1.28105615</td><td>0.04181181</td><td>-16.766942</td><td>4.258689e-63</td><td>1.514305e-59</td></tr>
	<tr><td>TRINITY_DN13736_c0_g1</td><td>NOC_DROME  </td><td>NOC  </td><td>  181.252900</td><td>-1.21565266</td><td>0.09412455</td><td> -6.753314</td><td>1.445052e-11</td><td>7.556349e-09</td></tr>
	<tr><td>TRINITY_DN13845_c0_g1</td><td>ABCA3_RAT  </td><td>ABCA3</td><td> 4634.329074</td><td> 0.25131029</td><td>0.03858003</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN1388_c1_g1 </td><td>NOC_DROME  </td><td>NOC  </td><td> 1933.800505</td><td>-1.22289281</td><td>0.05766902</td><td>-11.147975</td><td>7.325411e-29</td><td>1.302385e-25</td></tr>
	<tr><td>TRINITY_DN1464_c0_g1 </td><td>ABCA3_MOUSE</td><td>ABCA3</td><td> 4845.338166</td><td>-0.13338383</td><td>0.03048585</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN15745_c2_g1</td><td>SCPDL_BOVIN</td><td>SCPDL</td><td>   59.710049</td><td> 0.16274492</td><td>0.17771876</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN15745_c2_g2</td><td>SCPDL_HUMAN</td><td>SCPDL</td><td>   10.118123</td><td> 0.35916292</td><td>0.33638758</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN16169_c0_g1</td><td>NAAT1_DROME</td><td>NAAT1</td><td>  147.531120</td><td>-0.22811678</td><td>0.21024704</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN1718_c3_g1 </td><td>ABCA3_RAT  </td><td>ABCA3</td><td> 1175.622414</td><td> 0.26456046</td><td>0.05996947</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN1856_c1_g1 </td><td>CP4CU_BLAGE</td><td>CP4CU</td><td>  665.072599</td><td> 0.03358851</td><td>0.04623694</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN20906_c0_g1</td><td>CUD2_SCHGR </td><td>CUD2 </td><td>    1.891540</td><td> 0.18167683</td><td>1.26135038</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN21797_c0_g1</td><td>SCX_CHICK  </td><td>SCX  </td><td>    9.462207</td><td> 1.97133810</td><td>0.33857834</td><td>  4.109354</td><td>3.967682e-05</td><td>8.692614e-03</td></tr>
	<tr><td>TRINITY_DN2193_c2_g1 </td><td>USH_DROME  </td><td>USH  </td><td> 1513.849900</td><td>-1.02091399</td><td>0.03752914</td><td>-11.748577</td><td>7.181815e-32</td><td>1.418728e-28</td></tr>
	<tr><td>TRINITY_DN2223_c0_g1 </td><td>CHIT1_HUMAN</td><td>CHIT1</td><td>  111.094967</td><td> 1.54744884</td><td>0.14781876</td><td>  6.544831</td><td>5.956256e-11</td><td>3.025608e-08</td></tr>
	<tr><td>TRINITY_DN2457_c0_g1 </td><td>GGT1_RAT   </td><td>GGT1 </td><td> 4096.320992</td><td>-0.07389721</td><td>0.04235271</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN2597_c0_g1 </td><td>ABCA3_HUMAN</td><td>ABCA3</td><td>   29.873568</td><td>-0.46986741</td><td>0.23748681</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN27521_c0_g2</td><td>SCPDL_HUMAN</td><td>SCPDL</td><td>   42.155710</td><td>-0.13880490</td><td>0.19503741</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN27521_c1_g1</td><td>SCPDL_BOVIN</td><td>SCPDL</td><td>  376.170507</td><td>-0.23367861</td><td>0.14612829</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN27666_c0_g1</td><td>SCPDL_HUMAN</td><td>SCPDL</td><td>   28.099500</td><td> 2.21407823</td><td>0.26492908</td><td>  6.167984</td><td>6.916625e-10</td><td>3.074267e-07</td></tr>
	<tr><td>TRINITY_DN30325_c0_g1</td><td>CHIT1_MOUSE</td><td>CHIT1</td><td>   12.910170</td><td>-0.45130779</td><td>0.61014466</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>TRINITY_DN3046_c2_g1 </td><td>KEN1_CULQU </td><td>KEN1 </td><td>  784.860077</td><td>-0.81376173</td><td>0.05994472</td><td> -3.899621</td><td>9.634321e-05</td><td>1.841813e-02</td></tr>
	<tr><td>TRINITY_DN305_c0_g1  </td><td>GROU_DROME </td><td>GROU </td><td> 2258.155234</td><td>-1.13713764</td><td>0.04275220</td><td>-13.031788</td><td>8.069384e-39</td><td>1.793320e-35</td></tr>
	<tr><td>TRINITY_DN305_c0_g2  </td><td>GROU_DROME </td><td>GROU </td><td> 1387.108679</td><td>-0.23655659</td><td>0.04847316</td><td>  0.000000</td><td>1.000000e+00</td><td>1.000000e+00</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>TRINITY_DN3059_c0_g2 </td><td>GGT1_HUMAN </td><td>GGT1 </td><td> 486.627680</td><td>-0.26060212</td><td>0.05616718</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN31810_c1_g1</td><td>LIPT_MOUSE </td><td>LIPT </td><td>  58.423914</td><td>-1.35924261</td><td>0.18613481</td><td> -4.186442</td><td> 2.833612e-05</td><td> 6.542700e-03</td></tr>
	<tr><td>TRINITY_DN3481_c0_g1 </td><td>ABCA3_RAT  </td><td>ABCA3</td><td>5676.181397</td><td> 0.20238567</td><td>0.03723447</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN3598_c0_g2 </td><td>GGT1_HUMAN </td><td>GGT1 </td><td>   4.974651</td><td> 2.46727853</td><td>0.48030586</td><td>  3.929326</td><td> 8.518414e-05</td><td> 1.701673e-02</td></tr>
	<tr><td>TRINITY_DN3792_c0_g2 </td><td>TOLL7_DROME</td><td>TOLL7</td><td>1373.814697</td><td>-1.02396424</td><td>0.04776658</td><td> -9.294454</td><td> 1.479643e-20</td><td> 1.644161e-17</td></tr>
	<tr><td>TRINITY_DN421_c5_g1  </td><td>SMAD6_MOUSE</td><td>SMAD6</td><td>1559.955083</td><td>-2.26535350</td><td>0.04798355</td><td>-35.123570</td><td>2.944430e-270</td><td>5.234902e-266</td></tr>
	<tr><td>TRINITY_DN427_c1_g1  </td><td>CUD2_SCHGR </td><td>CUD2 </td><td>3475.820867</td><td>-0.03526084</td><td>0.06786224</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN431_c6_g1  </td><td>KI26L_DROME</td><td>KI26L</td><td>3889.921597</td><td>-0.81655091</td><td>0.02903284</td><td> -8.147700</td><td> 3.709108e-16</td><td> 2.997465e-13</td></tr>
	<tr><td>TRINITY_DN4393_c0_g1 </td><td>TRAF4_MOUSE</td><td>TRAF4</td><td>3331.565264</td><td> 0.83995647</td><td>0.05485586</td><td>  4.738900</td><td> 2.148812e-06</td><td> 6.586211e-04</td></tr>
	<tr><td>TRINITY_DN4418_c0_g1 </td><td>REG5_DROME </td><td>REG5 </td><td> 234.593915</td><td>-2.40733818</td><td>0.10422039</td><td>-17.533403</td><td> 7.965599e-69</td><td> 4.720679e-65</td></tr>
	<tr><td>TRINITY_DN45461_c0_g2</td><td>ABCA3_HUMAN</td><td>ABCA3</td><td> 559.186352</td><td> 0.27719642</td><td>0.06167226</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN47503_c0_g1</td><td>ABCA3_RAT  </td><td>ABCA3</td><td>   5.048087</td><td> 2.43696550</td><td>0.47377350</td><td>  3.919522</td><td> 8.872490e-05</td><td> 1.752711e-02</td></tr>
	<tr><td>TRINITY_DN4903_c0_g1 </td><td>CUD2_SCHGR </td><td>CUD2 </td><td> 649.198090</td><td>-1.54647213</td><td>0.06362851</td><td>-15.189294</td><td> 4.163831e-52</td><td> 1.233813e-48</td></tr>
	<tr><td>TRINITY_DN52430_c0_g2</td><td>PPAF3_HOLDI</td><td>PPAF3</td><td>  57.279697</td><td>-0.01976401</td><td>0.18366065</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN52870_c0_g2</td><td>LIPT_MOUSE </td><td>LIPT </td><td> 497.694032</td><td>-1.35336724</td><td>0.07185651</td><td>-10.762660</td><td> 5.165705e-27</td><td> 8.349189e-24</td></tr>
	<tr><td>TRINITY_DN556_c1_g1  </td><td>RIM2_DROME </td><td>RIM2 </td><td>2340.675187</td><td>-0.94752949</td><td>0.03426598</td><td>-10.725783</td><td> 7.703194e-27</td><td> 1.141292e-23</td></tr>
	<tr><td>TRINITY_DN57977_c0_g1</td><td>ABCA3_HUMAN</td><td>ABCA3</td><td> 760.163327</td><td> 0.15891391</td><td>0.06006830</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN6187_c1_g1 </td><td>SAM11_MOUSE</td><td>SAM11</td><td> 376.092462</td><td>-0.92842728</td><td>0.07038838</td><td> -4.950068</td><td> 7.418739e-07</td><td> 2.536495e-04</td></tr>
	<tr><td>TRINITY_DN640_c5_g1  </td><td>ATG2B_MOUSE</td><td>ATG2B</td><td>4979.530169</td><td> 0.83221949</td><td>0.03142300</td><td>  8.026589</td><td> 1.002204e-15</td><td> 7.747033e-13</td></tr>
	<tr><td>TRINITY_DN6656_c0_g1 </td><td>NAAT1_DROWI</td><td>NAAT1</td><td>2810.048294</td><td> 1.02636060</td><td>0.04472126</td><td>  9.980947</td><td> 1.846947e-23</td><td> 2.525913e-20</td></tr>
	<tr><td>TRINITY_DN6747_c0_g1 </td><td>CP4CU_BLAGE</td><td>CP4CU</td><td>  15.077818</td><td> 1.92956431</td><td>0.31700706</td><td>  4.257206</td><td> 2.069977e-05</td><td> 4.842384e-03</td></tr>
	<tr><td>TRINITY_DN6921_c1_g1 </td><td>KI26L_DROME</td><td>KI26L</td><td>3411.845160</td><td>-0.76670654</td><td>0.04659151</td><td> -4.007308</td><td> 6.141470e-05</td><td> 1.284579e-02</td></tr>
	<tr><td>TRINITY_DN75_c0_g1   </td><td>ABCA3_HUMAN</td><td>ABCA3</td><td>  69.701847</td><td> 1.41717270</td><td>0.18096003</td><td>  4.626285</td><td> 3.722826e-06</td><td> 1.018279e-03</td></tr>
	<tr><td>TRINITY_DN758_c0_g1  </td><td>NAAT1_DROVI</td><td>NAAT1</td><td> 688.857473</td><td>-0.31135020</td><td>0.04862254</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN76899_c0_g1</td><td>KI26L_DROME</td><td>KI26L</td><td> 364.777351</td><td>-0.85492326</td><td>0.05218310</td><td> -5.268434</td><td> 1.375922e-07</td><td> 5.096359e-05</td></tr>
	<tr><td>TRINITY_DN77323_c0_g2</td><td>ABCA3_HUMAN</td><td>ABCA3</td><td>1732.627809</td><td> 0.20285509</td><td>0.05308622</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN8757_c0_g1 </td><td>ABCA3_MOUSE</td><td>ABCA3</td><td>6847.002622</td><td> 0.21876317</td><td>0.03471430</td><td>  0.000000</td><td> 1.000000e+00</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN9426_c0_g1 </td><td>USH_DROME  </td><td>USH  </td><td>  40.372401</td><td>-0.86583332</td><td>0.17953722</td><td> -1.592056</td><td> 1.113721e-01</td><td> 1.000000e+00</td></tr>
	<tr><td>TRINITY_DN955_c1_g1  </td><td>TOLL7_DROME</td><td>TOLL7</td><td> 745.982054</td><td>-0.95210558</td><td>0.06122813</td><td> -6.077363</td><td> 1.221749e-09</td><td> 5.171778e-07</td></tr>
	<tr><td>TRINITY_DN9925_c0_g1 </td><td>DSX_DROME  </td><td>DSX  </td><td>  39.233396</td><td>-1.77158621</td><td>0.18420961</td><td> -6.468643</td><td> 9.888689e-11</td><td> 4.883639e-08</td></tr>
</tbody>
</table>




```R
intrest_sig_norm <- norm_count %>%
        filter(transcriptid %in% intrest_sig$transcriptid)
```


```R
joined_df <- intrest_sig %>%
  left_join(intrest_sig_norm, by = "transcriptid")

joined_df <- joined_df %>%
  select(-gene, -log2FoldChange, -lfcSE, -stat, -pvalue	, -padj, -baseMean)
```


```R
joined_df
```


<table class="dataframe">
<caption>A data.frame: 61 × 50</caption>
<thead>
	<tr><th scope=col>transcriptid</th><th scope=col>gene_simple</th><th scope=col>C_1_A</th><th scope=col>C_1_B</th><th scope=col>C_1_C</th><th scope=col>C_1_D</th><th scope=col>C_1_E</th><th scope=col>C_1_F</th><th scope=col>C_8_A</th><th scope=col>C_8_B</th><th scope=col>...</th><th scope=col>V_24_C</th><th scope=col>V_24_D</th><th scope=col>V_24_E</th><th scope=col>V_24_F</th><th scope=col>V_96_A</th><th scope=col>V_96_B</th><th scope=col>V_96_C</th><th scope=col>V_96_D</th><th scope=col>V_96_E</th><th scope=col>V_96_F</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>...</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TRINITY_DN0_c4_g2    </td><td>ABCA3</td><td> 5566.318</td><td> 5581.261</td><td> 5590.258</td><td> 5570.123</td><td> 5577.681</td><td> 5576.116</td><td> 5584.356</td><td> 5583.464</td><td>...</td><td> 5576.589</td><td> 5575.2190</td><td> 5578.3310</td><td> 5572.85600</td><td> 5581.704</td><td> 5579.812</td><td> 5582.125</td><td> 5587.861</td><td> 5582.355</td><td> 5576.563</td></tr>
	<tr><td>TRINITY_DN10814_c4_g1</td><td>LIPT </td><td>  125.481</td><td>  137.601</td><td>  145.342</td><td>  128.534</td><td>  134.242</td><td>  133.142</td><td>  141.040</td><td>  139.777</td><td>...</td><td>  133.594</td><td>  132.4650</td><td>  134.9550</td><td>  130.34900</td><td>  136.720</td><td>  135.297</td><td>  137.621</td><td>  142.502</td><td>  137.437</td><td>  132.719</td></tr>
	<tr><td>TRINITY_DN1083_c3_g1 </td><td>NAAT1</td><td> 2364.318</td><td> 2379.261</td><td> 2388.258</td><td> 2368.123</td><td> 2375.681</td><td> 2374.116</td><td> 2382.356</td><td> 2381.464</td><td>...</td><td> 2374.589</td><td> 2373.2190</td><td> 2376.3310</td><td> 2370.85600</td><td> 2379.704</td><td> 2377.812</td><td> 2380.125</td><td> 2385.861</td><td> 2380.355</td><td> 2374.563</td></tr>
	<tr><td>TRINITY_DN11417_c0_g1</td><td>FAT  </td><td>12058.318</td><td>12073.261</td><td>12082.258</td><td>12062.123</td><td>12069.681</td><td>12068.116</td><td>12076.356</td><td>12075.464</td><td>...</td><td>12068.589</td><td>12067.2190</td><td>12070.3310</td><td>12064.85600</td><td>12073.704</td><td>12071.812</td><td>12074.125</td><td>12079.861</td><td>12074.355</td><td>12068.563</td></tr>
	<tr><td>TRINITY_DN1151_c0_g1 </td><td>PPAF3</td><td> 2389.694</td><td> 2406.624</td><td> 2412.602</td><td> 2396.642</td><td> 2391.587</td><td> 2397.614</td><td> 2402.334</td><td> 2411.884</td><td>...</td><td> 2399.861</td><td> 2399.0546</td><td> 2405.4011</td><td> 2400.73579</td><td> 2390.308</td><td> 2395.947</td><td> 2397.995</td><td> 2401.030</td><td> 2400.545</td><td> 2399.010</td></tr>
	<tr><td>TRINITY_DN1178_c0_g1 </td><td>CHIT1</td><td> 6094.438</td><td> 6158.518</td><td> 6056.383</td><td> 6073.085</td><td> 6066.819</td><td> 6086.720</td><td> 6179.900</td><td> 6222.427</td><td>...</td><td> 5907.823</td><td> 5994.6154</td><td> 5924.1817</td><td> 5947.60443</td><td> 6009.720</td><td> 6017.143</td><td> 6037.329</td><td> 5998.056</td><td> 5977.739</td><td> 6005.358</td></tr>
	<tr><td>TRINITY_DN119_c79_g1 </td><td>ABCA3</td><td>  641.341</td><td>  656.281</td><td>  665.279</td><td>  645.140</td><td>  652.703</td><td>  651.141</td><td>  659.388</td><td>  658.485</td><td>...</td><td>  651.611</td><td>  650.2360</td><td>  653.3620</td><td>  647.88000</td><td>  656.728</td><td>  654.835</td><td>  657.139</td><td>  662.886</td><td>  657.376</td><td>  651.579</td></tr>
	<tr><td>TRINITY_DN1244_c3_g1 </td><td>LIPT </td><td> 1855.318</td><td> 1870.261</td><td> 1879.258</td><td> 1859.123</td><td> 1866.681</td><td> 1865.116</td><td> 1873.356</td><td> 1872.464</td><td>...</td><td> 1865.589</td><td> 1864.2190</td><td> 1867.3310</td><td> 1861.85600</td><td> 1870.704</td><td> 1868.812</td><td> 1871.125</td><td> 1876.861</td><td> 1871.355</td><td> 1865.563</td></tr>
	<tr><td>TRINITY_DN13736_c0_g1</td><td>NOC  </td><td>  664.330</td><td>  679.276</td><td>  688.275</td><td>  668.134</td><td>  675.695</td><td>  674.132</td><td>  682.383</td><td>  681.479</td><td>...</td><td>  674.603</td><td>  673.2320</td><td>  676.3470</td><td>  670.87200</td><td>  679.724</td><td>  677.829</td><td>  680.134</td><td>  685.880</td><td>  680.369</td><td>  674.572</td></tr>
	<tr><td>TRINITY_DN13845_c0_g1</td><td>ABCA3</td><td>  201.233</td><td>  215.561</td><td>  224.254</td><td>  204.621</td><td>  211.906</td><td>  210.429</td><td>  218.663</td><td>  217.772</td><td>...</td><td>  211.076</td><td>  209.6450</td><td>  212.6470</td><td>  207.36000</td><td>  215.823</td><td>  214.106</td><td>  216.493</td><td>  221.922</td><td>  216.504</td><td>  210.901</td></tr>
	<tr><td>TRINITY_DN1388_c1_g1 </td><td>NOC  </td><td> 1689.318</td><td> 1704.261</td><td> 1713.258</td><td> 1693.123</td><td> 1700.681</td><td> 1699.116</td><td> 1707.356</td><td> 1706.464</td><td>...</td><td> 1699.589</td><td> 1698.2190</td><td> 1701.3310</td><td> 1695.85600</td><td> 1704.704</td><td> 1702.812</td><td> 1705.125</td><td> 1710.861</td><td> 1705.355</td><td> 1699.563</td></tr>
	<tr><td>TRINITY_DN1464_c0_g1 </td><td>ABCA3</td><td> 5319.318</td><td> 5334.261</td><td> 5343.258</td><td> 5323.123</td><td> 5330.681</td><td> 5329.116</td><td> 5337.356</td><td> 5336.464</td><td>...</td><td> 5329.589</td><td> 5328.2190</td><td> 5331.3310</td><td> 5325.85600</td><td> 5334.704</td><td> 5332.812</td><td> 5335.125</td><td> 5340.861</td><td> 5335.355</td><td> 5329.563</td></tr>
	<tr><td>TRINITY_DN15745_c2_g1</td><td>SCPDL</td><td>  373.727</td><td>  388.583</td><td>  397.522</td><td>  377.392</td><td>  384.961</td><td>  383.410</td><td>  391.653</td><td>  390.778</td><td>...</td><td>  383.948</td><td>  382.5120</td><td>  385.6330</td><td>  380.17800</td><td>  388.998</td><td>  387.157</td><td>  389.457</td><td>  395.141</td><td>  389.674</td><td>  383.876</td></tr>
	<tr><td>TRINITY_DN15745_c2_g2</td><td>SCPDL</td><td>  315.912</td><td>  330.710</td><td>  339.633</td><td>  319.508</td><td>  327.080</td><td>  325.514</td><td>  333.786</td><td>  332.914</td><td>...</td><td>  326.090</td><td>  324.6580</td><td>  327.7620</td><td>  322.31400</td><td>  331.127</td><td>  329.299</td><td>  331.641</td><td>  337.276</td><td>  331.804</td><td>  326.025</td></tr>
	<tr><td>TRINITY_DN16169_c0_g1</td><td>NAAT1</td><td> 2360.318</td><td> 2375.261</td><td> 2384.258</td><td> 2364.123</td><td> 2371.681</td><td> 2370.116</td><td> 2378.356</td><td> 2377.464</td><td>...</td><td> 2370.589</td><td> 2369.2190</td><td> 2372.3310</td><td> 2366.85600</td><td> 2375.704</td><td> 2373.812</td><td> 2376.125</td><td> 2381.861</td><td> 2376.355</td><td> 2370.563</td></tr>
	<tr><td>TRINITY_DN1718_c3_g1 </td><td>ABCA3</td><td>   63.595</td><td>   66.670</td><td>   69.697</td><td>   66.468</td><td>   64.796</td><td>   66.651</td><td>   72.962</td><td>   68.679</td><td>...</td><td>   64.095</td><td>   65.2280</td><td>   65.2860</td><td>   64.02600</td><td>   59.962</td><td>   59.172</td><td>   60.739</td><td>   63.341</td><td>   60.551</td><td>   60.436</td></tr>
	<tr><td>TRINITY_DN1856_c1_g1 </td><td>CP4CU</td><td> 3922.318</td><td> 3911.911</td><td> 3941.338</td><td> 3908.074</td><td> 3921.695</td><td> 3925.449</td><td> 3905.984</td><td> 3927.284</td><td>...</td><td> 3918.914</td><td> 3931.2190</td><td> 3934.3310</td><td> 3823.88578</td><td> 3813.073</td><td> 3843.097</td><td> 3832.907</td><td> 3884.220</td><td> 3878.928</td><td> 3844.748</td></tr>
	<tr><td>TRINITY_DN20906_c0_g1</td><td>CUD2 </td><td>  132.335</td><td>  132.335</td><td>  132.335</td><td>  132.335</td><td>  132.335</td><td>  132.335</td><td>  141.040</td><td>  132.335</td><td>...</td><td>   74.420</td><td>  121.9974</td><td>  194.0288</td><td>   93.94599</td><td>  132.335</td><td>   71.432</td><td>  132.335</td><td>   76.379</td><td>  132.335</td><td>  132.335</td></tr>
	<tr><td>TRINITY_DN21797_c0_g1</td><td>SCX  </td><td> 1362.318</td><td> 1377.261</td><td> 1386.258</td><td> 1366.123</td><td> 1373.681</td><td> 1372.116</td><td> 1380.356</td><td> 1379.464</td><td>...</td><td> 1372.589</td><td> 1371.2190</td><td> 1374.3310</td><td> 1376.22069</td><td> 1377.704</td><td> 1376.221</td><td> 1378.125</td><td> 1383.861</td><td> 1378.355</td><td> 1372.563</td></tr>
	<tr><td>TRINITY_DN2193_c2_g1 </td><td>USH  </td><td> 3331.318</td><td> 3346.261</td><td> 3355.258</td><td> 3335.123</td><td> 3342.681</td><td> 3341.116</td><td> 3349.356</td><td> 3348.464</td><td>...</td><td> 3341.589</td><td> 3340.2190</td><td> 3343.3310</td><td> 3337.85600</td><td> 3346.704</td><td> 3344.812</td><td> 3347.125</td><td> 3352.861</td><td> 3347.355</td><td> 3341.563</td></tr>
	<tr><td>TRINITY_DN2223_c0_g1 </td><td>CHIT1</td><td> 7024.318</td><td> 7039.261</td><td> 7048.258</td><td> 7028.123</td><td> 7035.681</td><td> 7034.116</td><td> 7042.356</td><td> 7041.464</td><td>...</td><td> 7034.589</td><td> 7033.2190</td><td> 7036.3310</td><td> 7030.85600</td><td> 7039.704</td><td> 7037.812</td><td> 7040.125</td><td> 7045.861</td><td> 7040.355</td><td> 7034.563</td></tr>
	<tr><td>TRINITY_DN2457_c0_g1 </td><td>GGT1 </td><td> 4069.318</td><td> 4084.261</td><td> 4093.258</td><td> 4073.123</td><td> 4080.681</td><td> 4079.116</td><td> 4087.356</td><td> 4086.464</td><td>...</td><td> 4079.589</td><td> 4078.2190</td><td> 4081.3310</td><td> 4075.85600</td><td> 4084.704</td><td> 4082.812</td><td> 4085.125</td><td> 4090.861</td><td> 4085.355</td><td> 4079.563</td></tr>
	<tr><td>TRINITY_DN2597_c0_g1 </td><td>ABCA3</td><td> 5408.318</td><td> 5423.261</td><td> 5432.258</td><td> 5412.123</td><td> 5419.681</td><td> 5418.116</td><td> 5426.356</td><td> 5425.464</td><td>...</td><td> 5418.589</td><td> 5417.2190</td><td> 5420.3310</td><td> 5414.85600</td><td> 5423.704</td><td> 5421.812</td><td> 5424.125</td><td> 5429.861</td><td> 5424.355</td><td> 5418.563</td></tr>
	<tr><td>TRINITY_DN27521_c0_g2</td><td>SCPDL</td><td>  116.421</td><td>  127.955</td><td>  135.361</td><td>  119.479</td><td>  124.666</td><td>  123.689</td><td>  131.501</td><td>  130.141</td><td>...</td><td>  124.004</td><td>  123.0350</td><td>  125.3430</td><td>  120.98200</td><td>  126.727</td><td>  125.329</td><td>  127.653</td><td>  132.355</td><td>  127.408</td><td>  122.985</td></tr>
	<tr><td>TRINITY_DN27521_c1_g1</td><td>SCPDL</td><td>  463.558</td><td>  478.443</td><td>  487.409</td><td>  467.298</td><td>  474.849</td><td>  473.292</td><td>  481.554</td><td>  480.656</td><td>...</td><td>  473.788</td><td>  472.3890</td><td>  475.5120</td><td>  470.05800</td><td>  478.867</td><td>  476.994</td><td>  479.302</td><td>  485.022</td><td>  479.533</td><td>  473.735</td></tr>
	<tr><td>TRINITY_DN27666_c0_g1</td><td>SCPDL</td><td> 2159.318</td><td> 2174.261</td><td> 2183.258</td><td> 2163.123</td><td> 2170.681</td><td> 2169.116</td><td> 2177.356</td><td> 2176.464</td><td>...</td><td> 2169.589</td><td> 2168.2190</td><td> 2171.3310</td><td> 2165.85600</td><td> 2174.704</td><td> 2172.812</td><td> 2175.125</td><td> 2180.861</td><td> 2175.355</td><td> 2169.563</td></tr>
	<tr><td>TRINITY_DN30325_c0_g1</td><td>CHIT1</td><td> 1233.318</td><td> 1246.810</td><td> 1257.258</td><td> 1237.123</td><td> 1244.681</td><td> 1243.116</td><td> 1246.810</td><td> 1246.810</td><td>...</td><td> 1243.589</td><td> 1242.2190</td><td> 1245.3310</td><td> 1239.85600</td><td> 1248.704</td><td> 1246.812</td><td> 1249.125</td><td> 1254.861</td><td> 1249.355</td><td> 1243.563</td></tr>
	<tr><td>TRINITY_DN3046_c2_g1 </td><td>KEN1 </td><td> 2061.318</td><td> 2076.261</td><td> 2085.258</td><td> 2065.123</td><td> 2072.681</td><td> 2071.116</td><td> 2079.356</td><td> 2078.464</td><td>...</td><td> 2071.589</td><td> 2070.2190</td><td> 2073.3310</td><td> 2067.85600</td><td> 2076.704</td><td> 2074.812</td><td> 2077.125</td><td> 2082.861</td><td> 2077.355</td><td> 2071.563</td></tr>
	<tr><td>TRINITY_DN305_c0_g1  </td><td>GROU </td><td> 3153.318</td><td> 3168.261</td><td> 3177.258</td><td> 3157.123</td><td> 3164.681</td><td> 3163.116</td><td> 3171.356</td><td> 3170.464</td><td>...</td><td> 3163.589</td><td> 3162.2190</td><td> 3165.3310</td><td> 3159.85600</td><td> 3168.704</td><td> 3166.812</td><td> 3169.125</td><td> 3174.861</td><td> 3169.355</td><td> 3163.563</td></tr>
	<tr><td>TRINITY_DN305_c0_g2  </td><td>GROU </td><td> 3856.318</td><td> 3871.261</td><td> 3880.258</td><td> 3860.123</td><td> 3867.681</td><td> 3866.116</td><td> 3874.356</td><td> 3873.464</td><td>...</td><td> 3866.589</td><td> 3865.2190</td><td> 3868.3310</td><td> 3862.85600</td><td> 3871.704</td><td> 3869.812</td><td> 3872.125</td><td> 3877.861</td><td> 3872.355</td><td> 3866.563</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td></td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>TRINITY_DN3059_c0_g2 </td><td>GGT1 </td><td>1958.318</td><td>1973.261</td><td>1982.258</td><td>1962.123</td><td>1969.681</td><td>1968.1160</td><td>1976.356</td><td>1975.464</td><td>...</td><td>1968.589</td><td>1967.219</td><td>1970.331</td><td>1964.856</td><td>1973.704</td><td>1971.812</td><td>1974.125</td><td>1979.861</td><td>1974.355</td><td>1968.563</td></tr>
	<tr><td>TRINITY_DN31810_c1_g1</td><td>LIPT </td><td> 328.854</td><td> 343.674</td><td> 352.606</td><td> 332.467</td><td> 340.045</td><td> 338.4810</td><td> 346.743</td><td> 345.866</td><td>...</td><td> 339.051</td><td> 337.620</td><td> 340.726</td><td> 335.273</td><td> 344.088</td><td> 342.251</td><td> 344.582</td><td> 350.241</td><td> 344.770</td><td> 338.977</td></tr>
	<tr><td>TRINITY_DN3481_c0_g1 </td><td>ABCA3</td><td>4975.503</td><td>5142.407</td><td>5211.577</td><td>5107.731</td><td>5221.522</td><td>4946.1931</td><td>5168.657</td><td>5168.553</td><td>...</td><td>5105.443</td><td>5044.871</td><td>5137.481</td><td>5148.959</td><td>5111.490</td><td>5116.414</td><td>5225.433</td><td>5157.872</td><td>5083.602</td><td>5102.811</td></tr>
	<tr><td>TRINITY_DN3598_c0_g2 </td><td>GGT1 </td><td>2385.318</td><td>2400.261</td><td>2409.258</td><td>2389.123</td><td>2396.681</td><td>2399.0079</td><td>2403.356</td><td>2402.464</td><td>...</td><td>2395.589</td><td>2394.219</td><td>2399.008</td><td>2391.856</td><td>2399.008</td><td>2398.812</td><td>2401.125</td><td>2406.861</td><td>2401.355</td><td>2395.563</td></tr>
	<tr><td>TRINITY_DN3792_c0_g2 </td><td>TOLL7</td><td>8028.318</td><td>8043.261</td><td>8052.258</td><td>8032.123</td><td>8039.681</td><td>8038.1160</td><td>8046.356</td><td>8045.464</td><td>...</td><td>8038.589</td><td>8037.219</td><td>8040.331</td><td>8034.856</td><td>8043.704</td><td>8041.812</td><td>8044.125</td><td>8049.861</td><td>8044.355</td><td>8038.563</td></tr>
	<tr><td>TRINITY_DN421_c5_g1  </td><td>SMAD6</td><td>1153.318</td><td>1168.261</td><td>1177.258</td><td>1157.123</td><td>1164.681</td><td>1163.1160</td><td>1171.356</td><td>1170.464</td><td>...</td><td>1163.589</td><td>1162.219</td><td>1165.331</td><td>1159.856</td><td>1168.704</td><td>1166.812</td><td>1169.125</td><td>1174.861</td><td>1169.355</td><td>1163.563</td></tr>
	<tr><td>TRINITY_DN427_c1_g1  </td><td>CUD2 </td><td>1244.900</td><td>1312.243</td><td>1243.303</td><td>1305.536</td><td>1248.773</td><td>1236.1338</td><td>1076.778</td><td>1086.030</td><td>...</td><td>1376.271</td><td>1357.114</td><td>1357.201</td><td>1349.174</td><td>1445.840</td><td>1424.979</td><td>1426.058</td><td>1454.380</td><td>1433.611</td><td>1400.083</td></tr>
	<tr><td>TRINITY_DN431_c6_g1  </td><td>KI26L</td><td>1840.318</td><td>1855.261</td><td>1864.258</td><td>1844.123</td><td>1851.681</td><td>1850.1160</td><td>1858.356</td><td>1857.464</td><td>...</td><td>1850.589</td><td>1849.219</td><td>1852.331</td><td>1846.856</td><td>1855.704</td><td>1853.812</td><td>1856.125</td><td>1861.861</td><td>1856.355</td><td>1850.563</td></tr>
	<tr><td>TRINITY_DN4393_c0_g1 </td><td>TRAF4</td><td>1036.131</td><td>1071.757</td><td>1236.270</td><td>1122.668</td><td>1052.563</td><td>1056.3885</td><td>1163.578</td><td>1106.928</td><td>...</td><td>1071.688</td><td>1088.275</td><td>1084.460</td><td>1068.220</td><td>1078.803</td><td>1101.516</td><td>1013.358</td><td>1193.076</td><td>1105.697</td><td>1078.502</td></tr>
	<tr><td>TRINITY_DN4418_c0_g1 </td><td>REG5 </td><td> 995.318</td><td>1010.261</td><td>1019.258</td><td> 999.123</td><td>1006.681</td><td>1005.1160</td><td>1013.356</td><td>1012.464</td><td>...</td><td>1005.589</td><td>1004.219</td><td>1007.331</td><td>1001.856</td><td>1010.704</td><td>1008.812</td><td>1011.125</td><td>1016.861</td><td>1011.355</td><td>1005.563</td></tr>
	<tr><td>TRINITY_DN45461_c0_g2</td><td>ABCA3</td><td>  67.506</td><td>  71.863</td><td>  75.594</td><td>  70.417</td><td>  69.800</td><td>  71.1570</td><td>  77.624</td><td>  73.843</td><td>...</td><td>  69.065</td><td>  69.850</td><td>  70.311</td><td>  68.499</td><td>  65.959</td><td>  65.133</td><td>  66.766</td><td>  69.720</td><td>  66.559</td><td>  65.783</td></tr>
	<tr><td>TRINITY_DN47503_c0_g1</td><td>ABCA3</td><td> 218.820</td><td> 233.323</td><td> 242.074</td><td> 222.248</td><td> 229.651</td><td> 231.2175</td><td> 236.379</td><td> 235.503</td><td>...</td><td> 228.778</td><td> 227.341</td><td> 230.371</td><td> 225.019</td><td> 233.625</td><td> 231.875</td><td> 234.277</td><td> 239.743</td><td> 234.303</td><td> 228.648</td></tr>
	<tr><td>TRINITY_DN4903_c0_g1 </td><td>CUD2 </td><td>1726.318</td><td>1741.261</td><td>1750.258</td><td>1730.123</td><td>1737.681</td><td>1736.1160</td><td>1744.356</td><td>1743.464</td><td>...</td><td>1736.589</td><td>1735.219</td><td>1738.331</td><td>1732.856</td><td>1741.704</td><td>1739.812</td><td>1742.125</td><td>1747.861</td><td>1742.355</td><td>1736.563</td></tr>
	<tr><td>TRINITY_DN52430_c0_g2</td><td>PPAF3</td><td>3356.318</td><td>3371.261</td><td>3380.258</td><td>3360.123</td><td>3367.681</td><td>3366.1160</td><td>3374.356</td><td>3373.464</td><td>...</td><td>3366.589</td><td>3365.219</td><td>3368.331</td><td>3362.856</td><td>3371.704</td><td>3369.812</td><td>3372.125</td><td>3377.861</td><td>3372.355</td><td>3366.563</td></tr>
	<tr><td>TRINITY_DN52870_c0_g2</td><td>LIPT </td><td> 102.662</td><td> 112.952</td><td> 119.752</td><td> 105.592</td><td> 109.840</td><td> 109.1540</td><td> 116.803</td><td> 115.087</td><td>...</td><td> 109.178</td><td> 108.468</td><td> 110.445</td><td> 106.518</td><td> 110.936</td><td> 109.686</td><td> 111.916</td><td> 116.346</td><td> 111.608</td><td> 107.808</td></tr>
	<tr><td>TRINITY_DN556_c1_g1  </td><td>RIM2 </td><td>1035.318</td><td>1050.261</td><td>1059.258</td><td>1039.123</td><td>1046.681</td><td>1045.1160</td><td>1053.356</td><td>1052.464</td><td>...</td><td>1045.589</td><td>1044.219</td><td>1047.331</td><td>1041.856</td><td>1050.704</td><td>1048.812</td><td>1051.125</td><td>1056.861</td><td>1051.355</td><td>1045.563</td></tr>
	<tr><td>TRINITY_DN57977_c0_g1</td><td>ABCA3</td><td> 307.960</td><td> 322.736</td><td> 331.655</td><td> 311.537</td><td> 319.104</td><td> 317.5440</td><td> 325.810</td><td> 324.934</td><td>...</td><td> 318.122</td><td> 316.690</td><td> 319.790</td><td> 314.341</td><td> 323.153</td><td> 321.327</td><td> 323.673</td><td> 329.304</td><td> 323.830</td><td> 318.053</td></tr>
	<tr><td>TRINITY_DN6187_c1_g1 </td><td>SAM11</td><td>2082.074</td><td>2092.944</td><td>2096.755</td><td>2075.904</td><td>2093.424</td><td>2051.8057</td><td>2096.735</td><td>2095.719</td><td>...</td><td>2061.643</td><td>2073.959</td><td>2074.217</td><td>2071.672</td><td>2094.441</td><td>2085.706</td><td>2079.110</td><td>2070.651</td><td>2076.517</td><td>2095.563</td></tr>
	<tr><td>TRINITY_DN640_c5_g1  </td><td>ATG2B</td><td>6799.318</td><td>6814.261</td><td>6823.258</td><td>6803.123</td><td>6810.681</td><td>6809.1160</td><td>6817.356</td><td>6816.464</td><td>...</td><td>6809.589</td><td>6808.219</td><td>6811.331</td><td>6805.856</td><td>6814.704</td><td>6812.812</td><td>6815.125</td><td>6820.861</td><td>6815.355</td><td>6809.563</td></tr>
	<tr><td>TRINITY_DN6656_c0_g1 </td><td>NAAT1</td><td>2205.318</td><td>2220.261</td><td>2229.258</td><td>2209.123</td><td>2216.681</td><td>2215.1160</td><td>2223.356</td><td>2222.464</td><td>...</td><td>2215.589</td><td>2214.219</td><td>2217.331</td><td>2211.856</td><td>2220.704</td><td>2218.812</td><td>2221.125</td><td>2226.861</td><td>2221.355</td><td>2215.563</td></tr>
	<tr><td>TRINITY_DN6747_c0_g1 </td><td>CP4CU</td><td>5879.318</td><td>5894.261</td><td>5903.258</td><td>5883.123</td><td>5890.681</td><td>5889.1160</td><td>5897.356</td><td>5896.464</td><td>...</td><td>5889.589</td><td>5888.219</td><td>5891.331</td><td>5885.856</td><td>5894.704</td><td>5892.812</td><td>5895.125</td><td>5900.861</td><td>5895.355</td><td>5889.563</td></tr>
	<tr><td>TRINITY_DN6921_c1_g1 </td><td>KI26L</td><td>1654.318</td><td>1669.261</td><td>1678.258</td><td>1658.123</td><td>1665.681</td><td>1664.1160</td><td>1672.356</td><td>1671.464</td><td>...</td><td>1664.589</td><td>1663.219</td><td>1666.331</td><td>1660.856</td><td>1669.704</td><td>1667.812</td><td>1670.125</td><td>1675.861</td><td>1670.355</td><td>1664.563</td></tr>
	<tr><td>TRINITY_DN75_c0_g1   </td><td>ABCA3</td><td>4492.559</td><td>4436.768</td><td>4397.051</td><td>4542.767</td><td>4427.119</td><td>4561.1042</td><td>4706.356</td><td>4203.848</td><td>...</td><td>4390.533</td><td>4423.969</td><td>4343.046</td><td>4458.913</td><td>4132.213</td><td>4535.012</td><td>4704.125</td><td>4419.980</td><td>4329.734</td><td>4502.876</td></tr>
	<tr><td>TRINITY_DN758_c0_g1  </td><td>NAAT1</td><td>2998.318</td><td>3013.261</td><td>3022.258</td><td>3002.123</td><td>3009.681</td><td>3008.1160</td><td>3016.356</td><td>3015.464</td><td>...</td><td>3008.589</td><td>3007.219</td><td>3010.331</td><td>3004.856</td><td>3013.704</td><td>3011.812</td><td>3014.125</td><td>3019.861</td><td>3014.355</td><td>3008.563</td></tr>
	<tr><td>TRINITY_DN76899_c0_g1</td><td>KI26L</td><td> 217.840</td><td> 232.336</td><td> 241.081</td><td> 221.265</td><td> 228.661</td><td> 227.1570</td><td> 235.392</td><td> 234.516</td><td>...</td><td> 227.793</td><td> 226.355</td><td> 229.385</td><td> 224.034</td><td> 232.635</td><td> 230.886</td><td> 233.287</td><td> 238.751</td><td> 233.313</td><td> 227.659</td></tr>
	<tr><td>TRINITY_DN77323_c0_g2</td><td>ABCA3</td><td> 454.570</td><td> 469.457</td><td> 478.422</td><td> 458.308</td><td> 465.858</td><td> 464.3030</td><td> 472.562</td><td> 471.664</td><td>...</td><td> 464.794</td><td> 463.394</td><td> 466.522</td><td> 461.064</td><td> 469.877</td><td> 468.006</td><td> 470.318</td><td> 476.030</td><td> 470.541</td><td> 464.749</td></tr>
	<tr><td>TRINITY_DN8757_c0_g1 </td><td>ABCA3</td><td> 367.738</td><td> 382.595</td><td> 391.531</td><td> 371.401</td><td> 378.971</td><td> 377.4150</td><td> 385.664</td><td> 384.783</td><td>...</td><td> 377.959</td><td> 376.523</td><td> 379.644</td><td> 374.191</td><td> 383.007</td><td> 381.164</td><td> 383.468</td><td> 389.154</td><td> 383.688</td><td> 377.884</td></tr>
	<tr><td>TRINITY_DN9426_c0_g1 </td><td>USH  </td><td>  73.821</td><td>  79.700</td><td>  84.196</td><td>  76.646</td><td>  77.337</td><td>  78.1020</td><td>  84.825</td><td>  81.679</td><td>...</td><td>  76.573</td><td>  76.919</td><td>  77.850</td><td>  75.427</td><td>  74.767</td><td>  73.911</td><td>  75.626</td><td>  78.994</td><td>  75.408</td><td>  73.763</td></tr>
	<tr><td>TRINITY_DN955_c1_g1  </td><td>TOLL7</td><td>2165.318</td><td>2180.261</td><td>2189.258</td><td>2169.123</td><td>2176.681</td><td>2175.1160</td><td>2183.356</td><td>2182.464</td><td>...</td><td>2175.589</td><td>2174.219</td><td>2177.331</td><td>2171.856</td><td>2180.704</td><td>2178.812</td><td>2181.125</td><td>2186.861</td><td>2181.355</td><td>2175.563</td></tr>
	<tr><td>TRINITY_DN9925_c0_g1 </td><td>DSX  </td><td>2185.318</td><td>2200.261</td><td>2209.258</td><td>2189.123</td><td>2196.681</td><td>2195.1160</td><td>2203.356</td><td>2202.464</td><td>...</td><td>2195.589</td><td>2194.219</td><td>2197.331</td><td>2191.856</td><td>2200.704</td><td>2198.812</td><td>2201.125</td><td>2206.861</td><td>2201.355</td><td>2195.563</td></tr>
</tbody>
</table>




```R
# Gathering the columns to have normalized counts to a single column
joined_df_melt <- joined_df %>%
  gather(colnames(joined_df)[3:50], key = "SampleName", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(joined_df_melt)
```


<table class="dataframe">
<caption>A data.frame: 2928 × 4</caption>
<thead>
	<tr><th scope=col>transcriptid</th><th scope=col>gene_simple</th><th scope=col>SampleName</th><th scope=col>normalized_counts</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TRINITY_DN0_c4_g2    </td><td>ABCA3</td><td>C_1_A</td><td> 5566.318</td></tr>
	<tr><td>TRINITY_DN10814_c4_g1</td><td>LIPT </td><td>C_1_A</td><td>  125.481</td></tr>
	<tr><td>TRINITY_DN1083_c3_g1 </td><td>NAAT1</td><td>C_1_A</td><td> 2364.318</td></tr>
	<tr><td>TRINITY_DN11417_c0_g1</td><td>FAT  </td><td>C_1_A</td><td>12058.318</td></tr>
	<tr><td>TRINITY_DN1151_c0_g1 </td><td>PPAF3</td><td>C_1_A</td><td> 2389.694</td></tr>
	<tr><td>TRINITY_DN1178_c0_g1 </td><td>CHIT1</td><td>C_1_A</td><td> 6094.438</td></tr>
	<tr><td>TRINITY_DN119_c79_g1 </td><td>ABCA3</td><td>C_1_A</td><td>  641.341</td></tr>
	<tr><td>TRINITY_DN1244_c3_g1 </td><td>LIPT </td><td>C_1_A</td><td> 1855.318</td></tr>
	<tr><td>TRINITY_DN13736_c0_g1</td><td>NOC  </td><td>C_1_A</td><td>  664.330</td></tr>
	<tr><td>TRINITY_DN13845_c0_g1</td><td>ABCA3</td><td>C_1_A</td><td>  201.233</td></tr>
	<tr><td>TRINITY_DN1388_c1_g1 </td><td>NOC  </td><td>C_1_A</td><td> 1689.318</td></tr>
	<tr><td>TRINITY_DN1464_c0_g1 </td><td>ABCA3</td><td>C_1_A</td><td> 5319.318</td></tr>
	<tr><td>TRINITY_DN15745_c2_g1</td><td>SCPDL</td><td>C_1_A</td><td>  373.727</td></tr>
	<tr><td>TRINITY_DN15745_c2_g2</td><td>SCPDL</td><td>C_1_A</td><td>  315.912</td></tr>
	<tr><td>TRINITY_DN16169_c0_g1</td><td>NAAT1</td><td>C_1_A</td><td> 2360.318</td></tr>
	<tr><td>TRINITY_DN1718_c3_g1 </td><td>ABCA3</td><td>C_1_A</td><td>   63.595</td></tr>
	<tr><td>TRINITY_DN1856_c1_g1 </td><td>CP4CU</td><td>C_1_A</td><td> 3922.318</td></tr>
	<tr><td>TRINITY_DN20906_c0_g1</td><td>CUD2 </td><td>C_1_A</td><td>  132.335</td></tr>
	<tr><td>TRINITY_DN21797_c0_g1</td><td>SCX  </td><td>C_1_A</td><td> 1362.318</td></tr>
	<tr><td>TRINITY_DN2193_c2_g1 </td><td>USH  </td><td>C_1_A</td><td> 3331.318</td></tr>
	<tr><td>TRINITY_DN2223_c0_g1 </td><td>CHIT1</td><td>C_1_A</td><td> 7024.318</td></tr>
	<tr><td>TRINITY_DN2457_c0_g1 </td><td>GGT1 </td><td>C_1_A</td><td> 4069.318</td></tr>
	<tr><td>TRINITY_DN2597_c0_g1 </td><td>ABCA3</td><td>C_1_A</td><td> 5408.318</td></tr>
	<tr><td>TRINITY_DN27521_c0_g2</td><td>SCPDL</td><td>C_1_A</td><td>  116.421</td></tr>
	<tr><td>TRINITY_DN27521_c1_g1</td><td>SCPDL</td><td>C_1_A</td><td>  463.558</td></tr>
	<tr><td>TRINITY_DN27666_c0_g1</td><td>SCPDL</td><td>C_1_A</td><td> 2159.318</td></tr>
	<tr><td>TRINITY_DN30325_c0_g1</td><td>CHIT1</td><td>C_1_A</td><td> 1233.318</td></tr>
	<tr><td>TRINITY_DN3046_c2_g1 </td><td>KEN1 </td><td>C_1_A</td><td> 2061.318</td></tr>
	<tr><td>TRINITY_DN305_c0_g1  </td><td>GROU </td><td>C_1_A</td><td> 3153.318</td></tr>
	<tr><td>TRINITY_DN305_c0_g2  </td><td>GROU </td><td>C_1_A</td><td> 3856.318</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>TRINITY_DN3059_c0_g2 </td><td>GGT1 </td><td>V_96_F</td><td>1968.563</td></tr>
	<tr><td>TRINITY_DN31810_c1_g1</td><td>LIPT </td><td>V_96_F</td><td> 338.977</td></tr>
	<tr><td>TRINITY_DN3481_c0_g1 </td><td>ABCA3</td><td>V_96_F</td><td>5102.811</td></tr>
	<tr><td>TRINITY_DN3598_c0_g2 </td><td>GGT1 </td><td>V_96_F</td><td>2395.563</td></tr>
	<tr><td>TRINITY_DN3792_c0_g2 </td><td>TOLL7</td><td>V_96_F</td><td>8038.563</td></tr>
	<tr><td>TRINITY_DN421_c5_g1  </td><td>SMAD6</td><td>V_96_F</td><td>1163.563</td></tr>
	<tr><td>TRINITY_DN427_c1_g1  </td><td>CUD2 </td><td>V_96_F</td><td>1400.083</td></tr>
	<tr><td>TRINITY_DN431_c6_g1  </td><td>KI26L</td><td>V_96_F</td><td>1850.563</td></tr>
	<tr><td>TRINITY_DN4393_c0_g1 </td><td>TRAF4</td><td>V_96_F</td><td>1078.502</td></tr>
	<tr><td>TRINITY_DN4418_c0_g1 </td><td>REG5 </td><td>V_96_F</td><td>1005.563</td></tr>
	<tr><td>TRINITY_DN45461_c0_g2</td><td>ABCA3</td><td>V_96_F</td><td>  65.783</td></tr>
	<tr><td>TRINITY_DN47503_c0_g1</td><td>ABCA3</td><td>V_96_F</td><td> 228.648</td></tr>
	<tr><td>TRINITY_DN4903_c0_g1 </td><td>CUD2 </td><td>V_96_F</td><td>1736.563</td></tr>
	<tr><td>TRINITY_DN52430_c0_g2</td><td>PPAF3</td><td>V_96_F</td><td>3366.563</td></tr>
	<tr><td>TRINITY_DN52870_c0_g2</td><td>LIPT </td><td>V_96_F</td><td> 107.808</td></tr>
	<tr><td>TRINITY_DN556_c1_g1  </td><td>RIM2 </td><td>V_96_F</td><td>1045.563</td></tr>
	<tr><td>TRINITY_DN57977_c0_g1</td><td>ABCA3</td><td>V_96_F</td><td> 318.053</td></tr>
	<tr><td>TRINITY_DN6187_c1_g1 </td><td>SAM11</td><td>V_96_F</td><td>2095.563</td></tr>
	<tr><td>TRINITY_DN640_c5_g1  </td><td>ATG2B</td><td>V_96_F</td><td>6809.563</td></tr>
	<tr><td>TRINITY_DN6656_c0_g1 </td><td>NAAT1</td><td>V_96_F</td><td>2215.563</td></tr>
	<tr><td>TRINITY_DN6747_c0_g1 </td><td>CP4CU</td><td>V_96_F</td><td>5889.563</td></tr>
	<tr><td>TRINITY_DN6921_c1_g1 </td><td>KI26L</td><td>V_96_F</td><td>1664.563</td></tr>
	<tr><td>TRINITY_DN75_c0_g1   </td><td>ABCA3</td><td>V_96_F</td><td>4502.876</td></tr>
	<tr><td>TRINITY_DN758_c0_g1  </td><td>NAAT1</td><td>V_96_F</td><td>3008.563</td></tr>
	<tr><td>TRINITY_DN76899_c0_g1</td><td>KI26L</td><td>V_96_F</td><td> 227.659</td></tr>
	<tr><td>TRINITY_DN77323_c0_g2</td><td>ABCA3</td><td>V_96_F</td><td> 464.749</td></tr>
	<tr><td>TRINITY_DN8757_c0_g1 </td><td>ABCA3</td><td>V_96_F</td><td> 377.884</td></tr>
	<tr><td>TRINITY_DN9426_c0_g1 </td><td>USH  </td><td>V_96_F</td><td>  73.763</td></tr>
	<tr><td>TRINITY_DN955_c1_g1  </td><td>TOLL7</td><td>V_96_F</td><td>2175.563</td></tr>
	<tr><td>TRINITY_DN9925_c0_g1 </td><td>DSX  </td><td>V_96_F</td><td>2195.563</td></tr>
</tbody>
</table>




```R
##Importing metadata will be a slighly individual process.
##The important things to know is to 
##1. we want all cols to be as a factor and not Intergers/flots
##2. We need to make sure we have a column to build a basic model from, i.e. a column that contains the indivial groups
##3. Col 1 needs to be the same was what you ouput salmon folder name is

samples <- read.table(file = "metadata_numerictime.tsv.txt", header = TRUE) ## First we will import our metadata file
samples <- as.data.frame(unclass(samples),stringsAsFactors=TRUE) ## When importing the table we set all string colums to factors
samples$Time <- as.factor(samples$Time) ## If you have a colum that is numbers, you must manuall set it to a factor

# For this experiment we had one treatment factor across time, so we want a column that groups treatment and time
samples$Trt_Time <- paste(samples$Treatment, samples$Time, sep="_") #combines treatment column and time column seperate by _
samples$Trt_Time <- as.factor(samples$Trt_Time) # Setting this new column as a factor
samples ##Inspect your metadata to make sure it looks correct
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
joined_df_melt <- inner_join(samples, joined_df_melt)
```

    [1m[22mJoining with `by = join_by(SampleName)`
    


```R
joined_df_melt
```


<table class="dataframe">
<caption>A data.frame: 2928 × 7</caption>
<thead>
	<tr><th scope=col>SampleName</th><th scope=col>Treatment</th><th scope=col>Time</th><th scope=col>Trt_Time</th><th scope=col>transcriptid</th><th scope=col>gene_simple</th><th scope=col>normalized_counts</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN0_c4_g2    </td><td>ABCA3</td><td> 5566.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN10814_c4_g1</td><td>LIPT </td><td>  125.481</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1083_c3_g1 </td><td>NAAT1</td><td> 2364.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN11417_c0_g1</td><td>FAT  </td><td>12058.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1151_c0_g1 </td><td>PPAF3</td><td> 2389.694</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1178_c0_g1 </td><td>CHIT1</td><td> 6094.438</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN119_c79_g1 </td><td>ABCA3</td><td>  641.341</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1244_c3_g1 </td><td>LIPT </td><td> 1855.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN13736_c0_g1</td><td>NOC  </td><td>  664.330</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN13845_c0_g1</td><td>ABCA3</td><td>  201.233</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1388_c1_g1 </td><td>NOC  </td><td> 1689.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1464_c0_g1 </td><td>ABCA3</td><td> 5319.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN15745_c2_g1</td><td>SCPDL</td><td>  373.727</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN15745_c2_g2</td><td>SCPDL</td><td>  315.912</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN16169_c0_g1</td><td>NAAT1</td><td> 2360.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1718_c3_g1 </td><td>ABCA3</td><td>   63.595</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1856_c1_g1 </td><td>CP4CU</td><td> 3922.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN20906_c0_g1</td><td>CUD2 </td><td>  132.335</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN21797_c0_g1</td><td>SCX  </td><td> 1362.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2193_c2_g1 </td><td>USH  </td><td> 3331.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2223_c0_g1 </td><td>CHIT1</td><td> 7024.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2457_c0_g1 </td><td>GGT1 </td><td> 4069.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2597_c0_g1 </td><td>ABCA3</td><td> 5408.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN27521_c0_g2</td><td>SCPDL</td><td>  116.421</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN27521_c1_g1</td><td>SCPDL</td><td>  463.558</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN27666_c0_g1</td><td>SCPDL</td><td> 2159.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN30325_c0_g1</td><td>CHIT1</td><td> 1233.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN3046_c2_g1 </td><td>KEN1 </td><td> 2061.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN305_c0_g1  </td><td>GROU </td><td> 3153.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN305_c0_g2  </td><td>GROU </td><td> 3856.318</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN3059_c0_g2 </td><td>GGT1 </td><td>1968.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN31810_c1_g1</td><td>LIPT </td><td> 338.977</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN3481_c0_g1 </td><td>ABCA3</td><td>5102.811</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN3598_c0_g2 </td><td>GGT1 </td><td>2395.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN3792_c0_g2 </td><td>TOLL7</td><td>8038.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN421_c5_g1  </td><td>SMAD6</td><td>1163.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN427_c1_g1  </td><td>CUD2 </td><td>1400.083</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN431_c6_g1  </td><td>KI26L</td><td>1850.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN4393_c0_g1 </td><td>TRAF4</td><td>1078.502</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN4418_c0_g1 </td><td>REG5 </td><td>1005.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN45461_c0_g2</td><td>ABCA3</td><td>  65.783</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN47503_c0_g1</td><td>ABCA3</td><td> 228.648</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN4903_c0_g1 </td><td>CUD2 </td><td>1736.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN52430_c0_g2</td><td>PPAF3</td><td>3366.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN52870_c0_g2</td><td>LIPT </td><td> 107.808</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN556_c1_g1  </td><td>RIM2 </td><td>1045.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN57977_c0_g1</td><td>ABCA3</td><td> 318.053</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN6187_c1_g1 </td><td>SAM11</td><td>2095.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN640_c5_g1  </td><td>ATG2B</td><td>6809.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN6656_c0_g1 </td><td>NAAT1</td><td>2215.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN6747_c0_g1 </td><td>CP4CU</td><td>5889.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN6921_c1_g1 </td><td>KI26L</td><td>1664.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN75_c0_g1   </td><td>ABCA3</td><td>4502.876</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN758_c0_g1  </td><td>NAAT1</td><td>3008.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN76899_c0_g1</td><td>KI26L</td><td> 227.659</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN77323_c0_g2</td><td>ABCA3</td><td> 464.749</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN8757_c0_g1 </td><td>ABCA3</td><td> 377.884</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN9426_c0_g1 </td><td>USH  </td><td>  73.763</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN955_c1_g1  </td><td>TOLL7</td><td>2175.563</td></tr>
	<tr><td>V_96_F</td><td>Infected</td><td>96</td><td>Infected_96</td><td>TRINITY_DN9925_c0_g1 </td><td>DSX  </td><td>2195.563</td></tr>
</tbody>
</table>




```R
## Timepoint 1
joined_df_melt_Tp1 <- joined_df_melt %>%
        filter(Time == 1)
```


```R
joined_df_melt_Tp1
```


<table class="dataframe">
<caption>A data.frame: 732 × 7</caption>
<thead>
	<tr><th scope=col>SampleName</th><th scope=col>Treatment</th><th scope=col>Time</th><th scope=col>Trt_Time</th><th scope=col>transcriptid</th><th scope=col>gene_simple</th><th scope=col>normalized_counts</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN0_c4_g2    </td><td>ABCA3</td><td> 5566.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN10814_c4_g1</td><td>LIPT </td><td>  125.481</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1083_c3_g1 </td><td>NAAT1</td><td> 2364.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN11417_c0_g1</td><td>FAT  </td><td>12058.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1151_c0_g1 </td><td>PPAF3</td><td> 2389.694</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1178_c0_g1 </td><td>CHIT1</td><td> 6094.438</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN119_c79_g1 </td><td>ABCA3</td><td>  641.341</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1244_c3_g1 </td><td>LIPT </td><td> 1855.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN13736_c0_g1</td><td>NOC  </td><td>  664.330</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN13845_c0_g1</td><td>ABCA3</td><td>  201.233</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1388_c1_g1 </td><td>NOC  </td><td> 1689.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1464_c0_g1 </td><td>ABCA3</td><td> 5319.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN15745_c2_g1</td><td>SCPDL</td><td>  373.727</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN15745_c2_g2</td><td>SCPDL</td><td>  315.912</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN16169_c0_g1</td><td>NAAT1</td><td> 2360.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1718_c3_g1 </td><td>ABCA3</td><td>   63.595</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN1856_c1_g1 </td><td>CP4CU</td><td> 3922.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN20906_c0_g1</td><td>CUD2 </td><td>  132.335</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN21797_c0_g1</td><td>SCX  </td><td> 1362.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2193_c2_g1 </td><td>USH  </td><td> 3331.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2223_c0_g1 </td><td>CHIT1</td><td> 7024.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2457_c0_g1 </td><td>GGT1 </td><td> 4069.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN2597_c0_g1 </td><td>ABCA3</td><td> 5408.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN27521_c0_g2</td><td>SCPDL</td><td>  116.421</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN27521_c1_g1</td><td>SCPDL</td><td>  463.558</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN27666_c0_g1</td><td>SCPDL</td><td> 2159.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN30325_c0_g1</td><td>CHIT1</td><td> 1233.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN3046_c2_g1 </td><td>KEN1 </td><td> 2061.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN305_c0_g1  </td><td>GROU </td><td> 3153.318</td></tr>
	<tr><td>C_1_A</td><td>Uninfected</td><td>1</td><td>Uninfected_1</td><td>TRINITY_DN305_c0_g2  </td><td>GROU </td><td> 3856.318</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN3059_c0_g2 </td><td>GGT1 </td><td>1964.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN31810_c1_g1</td><td>LIPT </td><td> 335.178</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN3481_c0_g1 </td><td>ABCA3</td><td>5087.891</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN3598_c0_g2 </td><td>GGT1 </td><td>2391.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN3792_c0_g2 </td><td>TOLL7</td><td>8034.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN421_c5_g1  </td><td>SMAD6</td><td>1159.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN427_c1_g1  </td><td>CUD2 </td><td>1277.117</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN431_c6_g1  </td><td>KI26L</td><td>1846.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN4393_c0_g1 </td><td>TRAF4</td><td>1080.558</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN4418_c0_g1 </td><td>REG5 </td><td>1001.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN45461_c0_g2</td><td>ABCA3</td><td>  72.467</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN47503_c0_g1</td><td>ABCA3</td><td> 225.125</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN4903_c0_g1 </td><td>CUD2 </td><td>1732.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN52430_c0_g2</td><td>PPAF3</td><td>3362.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN52870_c0_g2</td><td>LIPT </td><td> 108.426</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN556_c1_g1  </td><td>RIM2 </td><td>1041.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN57977_c0_g1</td><td>ABCA3</td><td> 314.268</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN6187_c1_g1 </td><td>SAM11</td><td>2086.543</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN640_c5_g1  </td><td>ATG2B</td><td>6805.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN6656_c0_g1 </td><td>NAAT1</td><td>2211.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN6747_c0_g1 </td><td>CP4CU</td><td>5885.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN6921_c1_g1 </td><td>KI26L</td><td>1660.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN75_c0_g1   </td><td>ABCA3</td><td>4310.404</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN758_c0_g1  </td><td>NAAT1</td><td>3004.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN76899_c0_g1</td><td>KI26L</td><td> 224.143</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN77323_c0_g2</td><td>ABCA3</td><td> 460.928</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN8757_c0_g1 </td><td>ABCA3</td><td> 374.065</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN9426_c0_g1 </td><td>USH  </td><td>  78.935</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN955_c1_g1  </td><td>TOLL7</td><td>2171.680</td></tr>
	<tr><td>V_1_F</td><td>Infected</td><td>1</td><td>Infected_1</td><td>TRINITY_DN9925_c0_g1 </td><td>DSX  </td><td>2191.680</td></tr>
</tbody>
</table>




```R
## plot using ggplot2
ggplot(joined_df_melt_Tp1) +
        geom_point(aes(x = gene_simple, y = normalized_counts, color = Treatment)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Genes of Interest") +
        theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(plot.title = element_text(hjust = 0.5))
```


    
![png](output_20_0.png)
    



```R
volcplot_table <- read.csv(file = "deseq2_results/results_timepoint_1_swissnames_volcplot.csv", header = TRUE)
volcplot_table
```


<table class="dataframe">
<caption>A data.frame: 17779 × 10</caption>
<thead>
	<tr><th scope=col>X</th><th scope=col>transcriptid</th><th scope=col>gene</th><th scope=col>gene_simple</th><th scope=col>baseMean</th><th scope=col>log2FoldChange</th><th scope=col>lfcSE</th><th scope=col>stat</th><th scope=col>pvalue</th><th scope=col>padj</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>10459</td><td>TRINITY_DN421_c5_g1  </td><td>SMAD6_MOUSE </td><td>SMAD6       </td><td> 1559.95508</td><td>-2.2653535</td><td>0.04798355</td><td>-35.123570</td><td>2.94e-270</td><td>5.23e-266</td></tr>
	<tr><td>10474</td><td>TRINITY_DN4216_c3_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  498.93872</td><td>-2.2324698</td><td>0.06708391</td><td>-24.632878</td><td>5.62e-134</td><td>4.99e-130</td></tr>
	<tr><td>10930</td><td>TRINITY_DN4418_c0_g1 </td><td>REG5_DROME  </td><td>REG5        </td><td>  234.59392</td><td>-2.4073382</td><td>0.10422039</td><td>-17.533403</td><td> 7.97e-69</td><td> 4.72e-65</td></tr>
	<tr><td>17061</td><td>TRINITY_DN9008_c0_g2 </td><td>Unidentified</td><td>Unidentified</td><td>  298.39385</td><td> 2.3192683</td><td>0.10128988</td><td> 17.171194</td><td> 4.36e-66</td><td> 1.94e-62</td></tr>
	<tr><td> 1335</td><td>TRINITY_DN1244_c3_g1 </td><td>LIPT_MOUSE  </td><td>LIPT        </td><td>21101.54591</td><td>-1.2810561</td><td>0.04181181</td><td>-16.766942</td><td> 4.26e-63</td><td> 1.51e-59</td></tr>
	<tr><td>11925</td><td>TRINITY_DN4903_c0_g1 </td><td>CUD2_SCHGR  </td><td>CUD2        </td><td>  649.19809</td><td>-1.5464721</td><td>0.06362851</td><td>-15.189294</td><td> 4.16e-52</td><td> 1.23e-48</td></tr>
	<tr><td>12377</td><td>TRINITY_DN5183_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  883.10019</td><td> 1.7471950</td><td>0.08573741</td><td> 13.613602</td><td> 3.32e-42</td><td> 8.44e-39</td></tr>
	<tr><td> 7785</td><td>TRINITY_DN305_c0_g1  </td><td>GROU_DROME  </td><td>GROU        </td><td> 2258.15523</td><td>-1.1371376</td><td>0.04275220</td><td>-13.031788</td><td> 8.07e-39</td><td> 1.79e-35</td></tr>
	<tr><td> 5171</td><td>TRINITY_DN2193_c2_g1 </td><td>USH_DROME   </td><td>USH         </td><td> 1513.84990</td><td>-1.0209140</td><td>0.03752914</td><td>-11.748577</td><td> 7.18e-32</td><td> 1.42e-28</td></tr>
	<tr><td> 1964</td><td>TRINITY_DN1388_c1_g1 </td><td>NOC_DROME   </td><td>NOC         </td><td> 1933.80050</td><td>-1.2228928</td><td>0.05766902</td><td>-11.147975</td><td> 7.33e-29</td><td> 1.30e-25</td></tr>
	<tr><td>12536</td><td>TRINITY_DN52870_c0_g2</td><td>LIPT_MOUSE  </td><td>LIPT        </td><td>  497.69403</td><td>-1.3533672</td><td>0.07185651</td><td>-10.762660</td><td> 5.17e-27</td><td> 8.35e-24</td></tr>
	<tr><td>12954</td><td>TRINITY_DN556_c1_g1  </td><td>RIM2_DROME  </td><td>RIM2        </td><td> 2340.67519</td><td>-0.9475295</td><td>0.03426598</td><td>-10.725783</td><td> 7.70e-27</td><td> 1.14e-23</td></tr>
	<tr><td>14626</td><td>TRINITY_DN6656_c0_g1 </td><td>NAAT1_DROWI </td><td>NAAT1       </td><td> 2810.04829</td><td> 1.0263606</td><td>0.04472127</td><td>  9.980947</td><td> 1.85e-23</td><td> 2.53e-20</td></tr>
	<tr><td>15033</td><td>TRINITY_DN697_c15_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  328.74891</td><td>-1.2981687</td><td>0.07491873</td><td> -9.585969</td><td> 9.16e-22</td><td> 1.16e-18</td></tr>
	<tr><td> 8946</td><td>TRINITY_DN3566_c1_g1 </td><td>Unidentified</td><td>Unidentified</td><td>   50.71419</td><td>-2.3797728</td><td>0.18925628</td><td> -9.509712</td><td> 1.91e-21</td><td> 2.27e-18</td></tr>
	<tr><td> 9485</td><td>TRINITY_DN3792_c0_g2 </td><td>TOLL7_DROME </td><td>TOLL7       </td><td> 1373.81470</td><td>-1.0239642</td><td>0.04776658</td><td> -9.294454</td><td> 1.48e-20</td><td> 1.64e-17</td></tr>
	<tr><td> 1709</td><td>TRINITY_DN13321_c0_g1</td><td>ESMC_DROME  </td><td>ESMC        </td><td> 2333.35443</td><td>-1.0156062</td><td>0.04723667</td><td> -9.221779</td><td> 2.92e-20</td><td> 3.06e-17</td></tr>
	<tr><td> 7364</td><td>TRINITY_DN2892_c0_g1 </td><td>HXK2_DROME  </td><td>HXK2        </td><td>  411.48552</td><td>-1.2544923</td><td>0.07520840</td><td> -8.968310</td><td> 3.01e-19</td><td> 2.97e-16</td></tr>
	<tr><td>10699</td><td>TRINITY_DN4323_c0_g1 </td><td>ELBOW_DROME </td><td>ELBOW       </td><td> 1844.77563</td><td>-0.9576512</td><td>0.04272106</td><td> -8.839932</td><td> 9.58e-19</td><td> 8.96e-16</td></tr>
	<tr><td>  487</td><td>TRINITY_DN10814_c4_g1</td><td>LIPT_BOVIN  </td><td>LIPT        </td><td>  286.11856</td><td>-1.3470136</td><td>0.08709560</td><td> -8.806571</td><td> 1.29e-18</td><td> 1.15e-15</td></tr>
	<tr><td>16684</td><td>TRINITY_DN856_c3_g1  </td><td>NKD_AEDAE   </td><td>NKD         </td><td> 1371.01889</td><td>-0.9048226</td><td>0.03776523</td><td> -8.601102</td><td> 7.90e-18</td><td> 6.68e-15</td></tr>
	<tr><td>10674</td><td>TRINITY_DN431_c6_g1  </td><td>KI26L_DROME </td><td>KI26L       </td><td> 3889.92160</td><td>-0.8165509</td><td>0.02903284</td><td> -8.147700</td><td> 3.71e-16</td><td> 3.00e-13</td></tr>
	<tr><td>14294</td><td>TRINITY_DN640_c5_g1  </td><td>ATG2B_MOUSE </td><td>ATG2B       </td><td> 4979.53017</td><td> 0.8322195</td><td>0.03142300</td><td>  8.026589</td><td> 1.00e-15</td><td> 7.75e-13</td></tr>
	<tr><td>10862</td><td>TRINITY_DN4392_c1_g1 </td><td>Unidentified</td><td>Unidentified</td><td>  230.25550</td><td> 1.5955171</td><td>0.12759027</td><td>  7.959205</td><td> 1.73e-15</td><td> 1.28e-12</td></tr>
	<tr><td> 2434</td><td>TRINITY_DN1486_c2_g2 </td><td>LRK1_CAEEL  </td><td>LRK1        </td><td>23950.33322</td><td>-0.8004437</td><td>0.02784276</td><td> -7.917454</td><td> 2.42e-15</td><td> 1.72e-12</td></tr>
	<tr><td> 2521</td><td>TRINITY_DN1509_c1_g1 </td><td>PTP10_DROME </td><td>PTP10       </td><td> 6288.83536</td><td> 0.8216252</td><td>0.03071754</td><td>  7.866033</td><td> 3.66e-15</td><td> 2.50e-12</td></tr>
	<tr><td> 4811</td><td>TRINITY_DN2076_c0_g1 </td><td>Unidentified</td><td>Unidentified</td><td> 3115.08851</td><td>-0.8408772</td><td>0.03390632</td><td> -7.694060</td><td> 1.43e-14</td><td> 9.39e-12</td></tr>
	<tr><td> 9760</td><td>TRINITY_DN3888_c0_g1 </td><td>FGFR1_DROME </td><td>FGFR1       </td><td> 1265.62527</td><td>-0.9286297</td><td>0.04566059</td><td> -7.635243</td><td> 2.25e-14</td><td> 1.43e-11</td></tr>
	<tr><td> 9721</td><td>TRINITY_DN3872_c0_g1 </td><td>PP1RA_MOUSE </td><td>PP1RA       </td><td> 5309.05720</td><td> 0.9481387</td><td>0.04927873</td><td>  7.470539</td><td> 7.99e-14</td><td> 4.90e-11</td></tr>
	<tr><td>  861</td><td>TRINITY_DN11539_c0_g1</td><td>Unidentified</td><td>Unidentified</td><td> 7109.85841</td><td> 0.9212399</td><td>0.04639469</td><td>  7.355150</td><td> 1.91e-13</td><td> 1.13e-10</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>17749</td><td>TRINITY_DN994_c5_g1  </td><td>Unidentified</td><td></td><td>   3.667495</td><td> 0.05007459</td><td>0.45635233</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17750</td><td>TRINITY_DN9941_c0_g1 </td><td>CBPA1_ANOGA </td><td></td><td>  11.475510</td><td> 2.11227073</td><td>0.60992114</td><td> 2.51224400</td><td>0.01199661</td><td>1</td></tr>
	<tr><td>17751</td><td>TRINITY_DN99413_c0_g1</td><td>Unidentified</td><td></td><td>  12.710134</td><td> 0.15901563</td><td>0.31886841</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17752</td><td>TRINITY_DN9942_c0_g1 </td><td>Unidentified</td><td></td><td>   7.525042</td><td>-0.62233573</td><td>0.45227485</td><td>-0.09360619</td><td>0.92542199</td><td>1</td></tr>
	<tr><td>17753</td><td>TRINITY_DN99429_c0_g1</td><td>Unidentified</td><td></td><td>  90.545316</td><td>-0.04628323</td><td>0.13140714</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17754</td><td>TRINITY_DN9943_c0_g1 </td><td>Unidentified</td><td></td><td>  23.639197</td><td> 0.55942322</td><td>0.20922604</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17755</td><td>TRINITY_DN9946_c0_g2 </td><td>Unidentified</td><td></td><td>  10.267645</td><td> 0.15027193</td><td>0.39868859</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17756</td><td>TRINITY_DN9949_c0_g1 </td><td>Unidentified</td><td></td><td>  12.921730</td><td>-0.36038790</td><td>0.45263255</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17757</td><td>TRINITY_DN995_c0_g1  </td><td>Unidentified</td><td></td><td> 204.814695</td><td>-0.19168843</td><td>0.09196291</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17758</td><td>TRINITY_DN995_c1_g1  </td><td>Unidentified</td><td></td><td>   5.118045</td><td> 1.36122006</td><td>0.67251967</td><td> 1.16163153</td><td>0.24538517</td><td>1</td></tr>
	<tr><td>17759</td><td>TRINITY_DN9952_c0_g1 </td><td>Unidentified</td><td></td><td> 285.904262</td><td> 0.12306739</td><td>0.15787161</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17760</td><td>TRINITY_DN9953_c0_g2 </td><td>Unidentified</td><td></td><td>  88.264371</td><td> 0.09916924</td><td>0.18008229</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17761</td><td>TRINITY_DN996_c1_g1  </td><td>CUD4_LOCMI  </td><td></td><td>  40.733987</td><td> 0.17545432</td><td>0.47467798</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17762</td><td>TRINITY_DN9962_c0_g3 </td><td>Unidentified</td><td></td><td>  13.595194</td><td> 0.48665882</td><td>0.32965248</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17763</td><td>TRINITY_DN997_c0_g3  </td><td>Unidentified</td><td></td><td>  21.848769</td><td>-0.48491941</td><td>0.28825167</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17764</td><td>TRINITY_DN99759_c0_g1</td><td>Unidentified</td><td></td><td>  92.053494</td><td> 0.20530985</td><td>0.11842520</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17765</td><td>TRINITY_DN99790_c0_g1</td><td>SYTC_HUMAN  </td><td></td><td>   9.098529</td><td> 0.33553979</td><td>0.35024460</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17766</td><td>TRINITY_DN998_c0_g1  </td><td>Unidentified</td><td></td><td> 272.405776</td><td>-0.10482563</td><td>0.07448690</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17767</td><td>TRINITY_DN998_c1_g1  </td><td>AXN_DROME   </td><td></td><td>2507.792218</td><td> 0.07044966</td><td>0.03003691</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17768</td><td>TRINITY_DN9980_c0_g1 </td><td>PO210_DROME </td><td></td><td>  99.680791</td><td> 0.15200640</td><td>0.09706556</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17769</td><td>TRINITY_DN9984_c0_g1 </td><td>CASP6_CHICK </td><td></td><td> 189.712394</td><td> 0.15984120</td><td>0.09915854</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17770</td><td>TRINITY_DN9986_c1_g1 </td><td>Unidentified</td><td></td><td>   4.793399</td><td> 0.24424601</td><td>0.59766434</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17771</td><td>TRINITY_DN999_c0_g1  </td><td>Unidentified</td><td></td><td>  13.147643</td><td>-0.06295881</td><td>0.28812429</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17772</td><td>TRINITY_DN99904_c1_g1</td><td>Unidentified</td><td></td><td>   6.671768</td><td> 0.14528915</td><td>0.60054748</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17773</td><td>TRINITY_DN9991_c0_g1 </td><td>Unidentified</td><td></td><td>  57.855916</td><td> 0.12116419</td><td>0.15767282</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17774</td><td>TRINITY_DN9993_c0_g2 </td><td>Unidentified</td><td></td><td>  77.248675</td><td>-0.11538146</td><td>0.17962431</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17775</td><td>TRINITY_DN9994_c3_g1 </td><td>Unidentified</td><td></td><td>  43.071385</td><td>-0.07261690</td><td>0.19921284</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17776</td><td>TRINITY_DN9995_c0_g1 </td><td>Unidentified</td><td></td><td>   7.219157</td><td> 0.38877925</td><td>0.48482995</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17777</td><td>TRINITY_DN9999_c0_g1 </td><td>TDC1_CAEEL  </td><td></td><td>  37.655719</td><td> 0.11657106</td><td>0.18765409</td><td> 0.00000000</td><td>1.00000000</td><td>1</td></tr>
	<tr><td>17778</td><td>TRINITY_DN9999_c1_g1 </td><td>Unidentified</td><td></td><td>   9.670259</td><td>-0.76154938</td><td>0.38178856</td><td>-0.47552335</td><td>0.63441399</td><td>1</td></tr>
</tbody>
</table>




```R
if (!dir.exists("plots/Enhanced_volcano")){
  dir.create("plots/Enhanced_volcano")
}else{
  print("dir exists")
}



#, 'TRINITY_DN7897_c1_g1', 'TRINITY_DN5183_c0_g1','TRINITY_DN9008_c0_g2', 'TRINITY_DN3566_c1_g1'
# Set it globally:
options(repr.plot.width=12, repr.plot.height=8)

png("plots/Enhanced_volcano/Enhanced_volcano.png" , width = 4800, height = 3200, res=400)


EnhancedVolcano(result_table_1,
    lab = result_table_1$gene_simple,
    x = 'log2FoldChange',
    y = 'padj',
    legendLabels=c('Not sig.', ~Log[2]~ '(FC)','adjp', 'adjp &' ~Log[2]~ '(FC)'),
    title = " ",
    subtitle = " ",
    ylab = bquote(~-Log[10] ~ '(' ~ italic(adjp) ~ ')'),
    selectLab = c(genes_interest),
    pCutoff = 0.05,
    FCcutoff = 0.58,
    pointSize = 4.0,
    colAlpha = 4/5,
    labSize = 6.0,
    boxedLabels = FALSE,
    labCol = 'black',
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    colConnectors = 'black',
    col=c('grey', 'grey', 'grey', '#ff2c79ff')) +
    ggplot2::coord_cartesian(xlim=c(-3, 3), ylim=c(0,75)) +
    ggplot2::scale_x_continuous(breaks=seq(-3,3, 1))

dev.off()
```

    [1] "dir exists"
    

    [1m[22mCoordinate system already present. Adding new coordinate system, which will
    replace the existing one.
    [1m[22mScale for [32mx[39m is already present.
    Adding another scale for [32mx[39m, which will replace the existing scale.
    Warning message:
    "ggrepel: 29 unlabeled data points (too many overlaps). Consider increasing max.overlaps"
    


<strong>png:</strong> 2



```R

```
