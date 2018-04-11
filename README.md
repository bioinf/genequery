# GeneQuery database processing pipeline 

This is the pipeline that allows you to download many microarray GSE datasets, and turn them into a [GeneQuery](http://artyomovlab.wustl.edu/genequery/searcher/) database. Following steps are included: 

* Downloading and processing of GPL annotations;
* Searching for, and downloading GSE expression matrices;
* Preprocessing and quality control of expression matrices;
* Iterative clustering with WGCNA;
* Preparation of eigengene expression heatmaps.

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Maxim Artyomov Laboratory](https://artyomovlab.wustl.edu/site/), [Washington University in St Louis](https://wustl.edu/)

(c) 2014-2018, GPL v3 license

## GPL selection and processing 

Species-specific search and platform usage statistics allow us to identify the top microarray platforms used for every organism. We have downloaded and prepreprocessed 31 annotation file for human, 30 for mouse, and 23 for rat. These cover the vast majority of samples (over 90%). Processed files are available in this repository (see GPL directory). 

Should you require additional GPLs, there are three options available at NCBI: 
* *annot*-type annotaton file (not always available); 
* *soft*-type annotation file;
* *miniml* XML-formatted annotation file. 

Two latter ones also include information about all samples that are deposited into GEO under this GPL, and are thus quite bulky. First option (*annot*) is preferable. 

The resulting files are then converted into a *3col* file that contains 
* probe ID; 
* HGNC gene symbol; 
* Entrez gene ID (which could be used to pull all the relevant info from [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene). 

These *3col* files are then used for further processing of expression matrices. 

Scripts **annot_to_3col.sh** and **soft_to_3col.sh** are used to extract the necessary columns from the annotations. 

## Finding all of the GSE experiments for an organism

Given below is a sample query at [GEO datasets](https://www.ncbi.nlm.nih.gov/gds/) for human (for mouse/rat, *Organism* identifier needs to be appropriately replaced):

`("Homo Sapiens"[Organism]) AND expression profiling by array[DataSet Type] AND "gse"[Filter]`

(We used to include `AND (8[n_samples] : 200[n_samples])`, but *n_samples* does not seem to be supported by GEO anymore) 
The results are then saved using *Send to* button below the search bar; save it as a file of summaries in default order. Example search results are included in **GEO_search** folder in this repository. 

![alt text](img/GEO_search.png?raw=true "Sending all the found GSE results to a text file")

## GSE series matrix download

After you get the metafile downloaded (named *gds_list.txt* by default. You can use sample search results provided in **GEO_search**), simply run

`gse_download.sh gds_list.txt`

This will stard the download of files named **GSEXXX.series_matrix.gz** 

The overall volume of files is ~ 20Gb for human, ~ 10Gb for mouse, and ~2Gb for rat. Depending on the Internet connectivity, this would take anywhere from several hours to a day. 

## QC and preprocessing 

Most of the QC and preprocessing are done by **preprocess_matrix.pl** Perl script. It takes two arguments: 
* a GSE series matrix file; 
* directory with GPL *3col* files (assuming they are sorted into *hs*, *mm*, and *rt* subdirectories). 

The script would attempt to guess if the matrix is log2 transformed or not (see below). The script would stop processing the matrix (consider QC failed) if one of the following conditions is satisfied: 
* GPL listed in the matrix file is absent in the directory of *3col* files; 
* Number of samples is fewer than a pre-defined minimum (default = 8); 
* Number of samples is greated than a pre-defined maximum (default = 200). 
* Expression matrix is empty; 
* Expression matrix maximum and average values are outside of pre-defined ranges: 
  * ave 4-10, max 10-25 for log2-transformed matrix;
  * max 5000-100000 for non-log2-transformed matrix. 

If none of the above criteria are satisfied, the matrix statistics are printed to logs (STDERR), and the non-log2-transformed matrix is printed following the 3col annotation columns in the comma-separated (CSV) format. 

## WGCNA clustering 

Scripts for WGCNA processing and heatmap generation are located in genequery/R. You need to have R and WGCNA installed. 

WGCNA clustering of top N genes by expression uses **GSEXXX_GPLYYY_preprocessed.csv**, generated as mentioned above. All of the processing is done by function *wgcna_preprocessed*, that does the following: 

* collapses probes to genes using WGCNA built-in function collapseRows;
* takes top N genes (6000 by default) as ranked by max2 function - see max2 definition in *genequery_functions.R*; 
* generates Rsq of scale-free fit for powers 5-25, after which
  * if Rsq is 0.9 is achievable, it picks the lowest power at which Rsq is >= 0.9;
  * if only Rsq of 0.8 is achievable, it picks the lowest power at which Rsq is >= 0.8; 
  * if only Rsq of 0.7 is achievable, it picks the lowest power at which Rsq is >= 0.7;
  * if Rsq is < 0.7 for all powers in [5:25], clustering is considered failed.
* after this, the function picks mergeCutHeight: 0.15 for pow > 10, and 0.25 for pow <= 10;
* finally, the clustering is performed, and several files are saved.

We would need genes-to-clusters tables (for GMT file that later becomes GeneQuery database), and eigengene table (for the heatmap). 

## Heatmap generation 

Heatmap generation is done by function *make_svg_heatmaps* in *genequery_functions.R*. It uses eigengene matrix generated during the WGCNA clustering. 

As a result, it generates N figures (one for each generated WGCNA cluster) in which cluster number *i* is highlighted by frame. Figures are saved in SVG format. 

