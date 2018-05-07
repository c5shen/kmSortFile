# kmSortFile
k-means clustering sorting to heatmap display

## Usage
1. This is an add-on to the **KMeanClustering** module from [GenePattern](https://genepattern.broadinstitute.org/gp/pages/login.jsf).
2. It will take in all output .gct files from the module (for each clusters and for the sample as a whole), and return a sorted output file.
3. The columns are ordered from cluster #1 to cluster #k, as outputted by the **KMeanClustering** module. The rows are sorted out by the level of differently expressed genes in clusters (default output by the module will not sort the rows and they will be in alphabetical order).
4. More specifically to *3*, the level of differently expressed genes is tested out using *two-sample t test* (with assumption of same variance). For each cluster from left to right (#1 to #k in the output files), each gene expression of its samples are  compared to the rest samples, and a p-value is recorded. If this p-value exceeds a threshold (default setting is 0.001), then we flag this gene as "differently expressed" on this particular cluster. The test is conducted for all genes for the first cluster, and will move on to the next cluster, and so on. For each time there exists a gene that is differently expressed in two (or more) different clusters, the p-values are compared and the smallest one is selected.

## Command
### Developed using Python 3.6.5
```bash
python3 kmSortFile.py [integrated output file name (no extension)] [number of clusters] [p value for t test, default: 0.005] [whether only output differentially expressed genes, default: 1]
```

**Example (files provided in the repository)**
```bash
python3 kmSortFile.py all_aml_test.preprocessed_KMcluster_output 4 0.005 1
```

**Output:**
a file named "all_aml_test.preprocessed_KMcluster_output-sorted.gct"
