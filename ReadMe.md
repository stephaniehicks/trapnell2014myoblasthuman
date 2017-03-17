# trapnell2014myoblasthuman

This package contains a Bioconductor 
    ExpressionSet from Trapnell et al. (2014) paper 
    that performed a time-series experiment of bulk 
    (three bulk samples per time point) and 
    single cell RNA-Seq at four time points (hour 0, 24, 48, 72) in 
    differentiated primary human myoblasts ([PMID: 24658644](http://www.ncbi.nlm.nih.gov/pubmed/24658644)). 
    Metadata and pre-processed data (FPKM) were downloaded
    from Gene Expression Omnibus ([GSE52529](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52529)). 
    Batch information was extracted from headers 
    of FASTQ files downloaded from Sequence Read Archive. 
    
# Installation

The R-package **trapnell2014myoblasthuman** can be installed from Github using the R
package **devtools**.
```s
library(devtools)
install_github("stephaniehicks/trapnell2014myoblasthuman")
```
# Load data

The data is provided as a `ExpressionSet` object can be loaded 
by running the following code in R: 

```r
library(trapnell2014myoblasthuman)
data(trapnell2014myoblasthuman)

# Get the expression data
eset = exprs(trapnell2014myoblasthuman)

# Get the pheno data
pd = pData(trapnell2014myoblasthuman)
```

# Bug reports
Report bugs as issues on the [GitHub repository](https://github.com/stephaniehicks/trapnell2014myoblasthuman)

# Contributors

* [Will Townes](https://github.com/willtownes)
* [Stephanie Hicks](https://github.com/stephaniehicks)