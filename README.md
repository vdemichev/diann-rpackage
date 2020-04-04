### DiaNN R Package
A package containing a number of convenient functions for dealing with feature quantification data. The package has been developed primarily to support DIA/SWATH-MS data analysis with DIA-NN https://github.com/vdemichev/DiaNN and allow for MaxLFQ-based protein quantification (https://doi.org/10.1074/mcp.M113.031591) after manual precursor-level filtering and batch correction.

**DIA-NN: neural networks and interference correction   
enable deep proteome coverage in high throughput**  
Vadim Demichev, Christoph B. Messner, Spyros I. Vernardis, Kathryn S. Lilley, Markus Ralser  
https://www.nature.com/articles/s41592-019-0638-x

**Installation in R**:
```
install.packages("devtools")
library(devtools)
install_github("https://github.com/vdemichev/diann-rpackage")
library(diann)
```

**Examples**:   
```
library(diann)
```
Load a DIA-NN report (a small sample report is included in this repository)
```
df <- diann_load("diann_report.tsv")
```
Precursors x samples matrix filtered at 1% precursor and protein group FDR
```
precursors <- diann_matrix(df, pg.q = 0.01)
```
Peptides without modifications - taking the maximum of the respective precursor quantities
```
peptides <- diann_matrix(df, id.header="Stripped.Sequence", pg.q = 0.01)
```
Peptides without modifications - using the MaxLFQ algorithm
```
peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Stripped.Sequence", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
```
Genes identified and quantified using proteotypic peptides
```
unique.genes <- diann_matrix(df, id.header="Genes", quantity.header="Genes.MaxLFQ.Unique", proteotypic.only = T, pg.q = 0.01)
```
Protein group quantities using MaxLFQ algorithm
```
protein.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")
```
