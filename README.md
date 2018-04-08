# genequery processing pipeline 

This is the pipeline that allows you to download many GSE datasets, and turn them into a [GeneQuery](http://artyomovlab.wustl.edu/genequery/searcher/) database. Following steps are included: 

* Downloading and processing of GPL annotations;
* Searching for, and downloading preprocessed GSE expression matrices;
* Preprocessing and quality control of expression matrices;
* Iterative clustering with WGCNA;
* Preparation of eigengene expression heatmaps.

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Maxim Artyomov Laboratory](https://artyomovlab.wustl.edu/site/), [Washington University in St Louis](https://wustl.edu/)

(c) 2014-2018, GPL v3 license
GeneQuery DB preprocessing pipeline.

## GPL selection and processing 

You choose and process GPLs

## data download

Given below is a sample query at [GEO datasets](https://www.ncbi.nlm.nih.gov/gds/) for human:
`("Homo Sapiens"[Organism]) AND expression profiling by array[DataSet Type] AND "gse"[Filter]`

(We used to include `AND (8[n_samples] : 200[n_samples])`, but *n_samples* does not seem to be supported by GEO anymore) 
The results are then saved using *Send to* button below the search bar; save it as a file of summaries in default order. 

Example search results are included in **GEO_search** folder in this repository. 

## download

After you get the files downloaded (you can choose provided as well), simply run

`gse_download.sh gds_list.txt`

This will stard the download of files named **GSEXXX.series_matrix.gz** 

The overall volume of files is ~ 20Gb for human, ~ 10Gb for mouse, and ~2Gb for rat. 

## QC and preprocessing 

## WGCNA clustering 

## heatmap generation 

