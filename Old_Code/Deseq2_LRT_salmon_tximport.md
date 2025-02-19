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
```

    package 'RColorBrewer' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\edwardbird\AppData\Local\Temp\RtmpeqRMGd\downloaded_packages
    


```R
citation()
```


    
    To cite R in publications use:
    
      R Core Team (2022). R: A language and environment for statistical
      computing. R Foundation for Statistical Computing, Vienna, Austria.
      URL https://www.R-project.org/.
    
    A BibTeX entry for LaTeX users is
    
      @Manual{,
        title = {R: A Language and Environment for Statistical Computing},
        author = {{R Core Team}},
        organization = {R Foundation for Statistical Computing},
        address = {Vienna, Austria},
        year = {2022},
        url = {https://www.R-project.org/},
      }
    
    We have invested a lot of time and effort in creating R, please cite it
    when using it for data analysis. See also 'citation("pkgname")' for
    citing R packages.
    



```R
library(tximportData)
library(tximport)
library(DESeq2)
library(tidyverse)
library(DEGreport)
library(pheatmap)
library(RColorBrewer)
```


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
meta = as.data.frame.matrix(samples)
meta
```


<table class="dataframe">
<caption>A data.frame: 48 × 4</caption>
<thead>
	<tr><th></th><th scope=col>SampleName</th><th scope=col>Treatment</th><th scope=col>Time</th><th scope=col>Trt_Time</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>C_1_A </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><th scope=row>2</th><td>C_1_B </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><th scope=row>3</th><td>C_1_C </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><th scope=row>4</th><td>C_1_D </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><th scope=row>5</th><td>C_1_E </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><th scope=row>6</th><td>C_1_F </td><td>Uninfected</td><td>1 </td><td>Uninfected_1 </td></tr>
	<tr><th scope=row>7</th><td>C_8_A </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><th scope=row>8</th><td>C_8_B </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><th scope=row>9</th><td>C_8_C </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><th scope=row>10</th><td>C_8_D </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><th scope=row>11</th><td>C_8_E </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><th scope=row>12</th><td>C_8_F </td><td>Uninfected</td><td>8 </td><td>Uninfected_8 </td></tr>
	<tr><th scope=row>13</th><td>C_24_A</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><th scope=row>14</th><td>C_24_B</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><th scope=row>15</th><td>C_24_C</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><th scope=row>16</th><td>C_24_D</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><th scope=row>17</th><td>C_24_E</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><th scope=row>18</th><td>C_24_F</td><td>Uninfected</td><td>24</td><td>Uninfected_24</td></tr>
	<tr><th scope=row>19</th><td>C_96_A</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><th scope=row>20</th><td>C_96_B</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><th scope=row>21</th><td>C_96_C</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><th scope=row>22</th><td>C_96_D</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><th scope=row>23</th><td>C_96_E</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><th scope=row>24</th><td>C_96_F</td><td>Uninfected</td><td>96</td><td>Uninfected_96</td></tr>
	<tr><th scope=row>25</th><td>V_1_A </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><th scope=row>26</th><td>V_1_B </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><th scope=row>27</th><td>V_1_C </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><th scope=row>28</th><td>V_1_D </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><th scope=row>29</th><td>V_1_E </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><th scope=row>30</th><td>V_1_F </td><td>Infected  </td><td>1 </td><td>Infected_1   </td></tr>
	<tr><th scope=row>31</th><td>V_8_A </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><th scope=row>32</th><td>V_8_B </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><th scope=row>33</th><td>V_8_C </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><th scope=row>34</th><td>V_8_D </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><th scope=row>35</th><td>V_8_E </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><th scope=row>36</th><td>V_8_F </td><td>Infected  </td><td>8 </td><td>Infected_8   </td></tr>
	<tr><th scope=row>37</th><td>V_24_A</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><th scope=row>38</th><td>V_24_B</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><th scope=row>39</th><td>V_24_C</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><th scope=row>40</th><td>V_24_D</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><th scope=row>41</th><td>V_24_E</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><th scope=row>42</th><td>V_24_F</td><td>Infected  </td><td>24</td><td>Infected_24  </td></tr>
	<tr><th scope=row>43</th><td>V_96_A</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><th scope=row>44</th><td>V_96_B</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><th scope=row>45</th><td>V_96_C</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><th scope=row>46</th><td>V_96_D</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><th scope=row>47</th><td>V_96_E</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
	<tr><th scope=row>48</th><td>V_96_F</td><td>Infected  </td><td>96</td><td>Infected_96  </td></tr>
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
dds <- DESeqDataSetFromTximport(txi.salmon,
                                   colData = samples,
                                   design = ~ Treatment + Time)
```

    using counts and average transcript lengths from tximport
    
    


```R
dds
```


    class: DESeqDataSet 
    dim: 48247 48 
    metadata(1): version
    assays(2): counts avgTxLength
    rownames(48247): TRINITY_DN0_c0_g1 TRINITY_DN0_c0_g2 ...
      TRINITY_DN9999_c0_g1 TRINITY_DN9999_c1_g1
    rowData names(0):
    colnames(48): C_1_A C_1_B ... V_96_E V_96_F
    colData names(4): SampleName Treatment Time Trt_Time



```R
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
```


```R
dds$condition <- relevel(dds$Treatment, ref = "Uninfected")
```


```R
dds <- DESeq(dds, test="LRT", reduced = ~ Time)
res_LRT <- results(dds)
```

    estimating size factors
    
    using 'avgTxLength' from assays(dds), correcting for library size
    
    estimating dispersions
    
    gene-wise dispersion estimates
    
    mean-dispersion relationship
    
    final dispersion estimates
    
    fitting model and testing
    
    


```R
res_LRT
```


    log2 fold change (MLE): Time 96 vs 1 
    LRT p-value: '~ Treatment + Time' vs '~ Time' 
    DataFrame with 17779 rows and 6 columns
                           baseMean log2FoldChange     lfcSE        stat
                          <numeric>      <numeric> <numeric>   <numeric>
    TRINITY_DN0_c0_g1    4898.71210      -0.433740 0.0274586  4.68775865
    TRINITY_DN0_c0_g2       9.29955      -0.524085 0.2605267  0.00339354
    TRINITY_DN0_c1_g1     259.54213       1.932209 0.0668125 13.77454347
    TRINITY_DN0_c10_g1   1446.01769       0.225846 0.0248291 14.77996670
    TRINITY_DN0_c11_g1      6.57610       0.101634 0.2900479  0.00365167
    ...                         ...            ...       ...         ...
    TRINITY_DN9993_c0_g2   77.24867       0.753051  0.130352  0.00343556
    TRINITY_DN9994_c3_g1   43.07138       0.142159  0.141260  0.00351279
    TRINITY_DN9995_c0_g1    7.21916       1.522523  0.307472  0.00349240
    TRINITY_DN9999_c0_g1   37.65572      -0.739475  0.150070  0.27333190
    TRINITY_DN9999_c1_g1    9.67026       0.839482  0.254916  0.27376410
                              pvalue       padj
                           <numeric>  <numeric>
    TRINITY_DN0_c0_g1    0.030378249 0.15937756
    TRINITY_DN0_c0_g2    0.953546241         NA
    TRINITY_DN0_c1_g1    0.000206110 0.00703137
    TRINITY_DN0_c10_g1   0.000120812 0.00519938
    TRINITY_DN0_c11_g1   0.951813917         NA
    ...                          ...        ...
    TRINITY_DN9993_c0_g2    0.953260   0.979953
    TRINITY_DN9994_c3_g1    0.952738   0.979816
    TRINITY_DN9995_c0_g1    0.952875         NA
    TRINITY_DN9999_c0_g1    0.601105   0.812800
    TRINITY_DN9999_c1_g1    0.600818         NA



```R
# Subset the LRT results to return genes with padj < 0.05
sig_res_LRT <- res_LRT %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>% 
               as_tibble() %>% 
               filter(padj < 0.05)
 
# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
                pull(gene)
                
length(sigLRT_genes)
```


1330



```R
 # Transform counts for data visualization
 rld <- rlog(dds, blind=TRUE)
	
 # Plot PCA 
 plotPCA(rld, intgroup="condition")
	
 # Extract the rlog matrix from the object and compute pairwise correlation values
 rld_mat <- assay(rld)
 rld_cor <- cor(rld_mat)
```

    rlog() may take a few minutes with 30 or more samples,
    vst() is a much faster transformation
    
    


    
![png](output_14_1.png)
    



```R
 # Plot heatmap
df <- as.data.frame(colData(dds)[,c("Treatment","Time")])
pheatmap(rld_cor, annotation_col=df)
```


    
![png](output_15_0.png)
    



```R
# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[sig_res_LRT$gene, ]
```


```R
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = df, time = "Time", col=NULL)
```

    Working with 1330 genes.
    
    Working with 1330 genes after filtering: minc > 15
    
    [1m[22mJoining with `by = join_by(merge)`
    [1m[22mJoining with `by = join_by(merge)`
    


    
![png](output_17_1.png)
    



```R

```
