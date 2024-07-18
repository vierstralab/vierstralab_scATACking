# vierstralab_scATACking
In this repository you can find python scrips for scATAC-seq data analysis with excamples of data and graphs. Also you may find the conda enviroment with all of the dependencies needed. Before working with thethese scripts, you should do the following steps:

1. Exclude ENCODE blacklist regions from each of the fragment files of a sample and intersecting the resulting files with the GRCh38 version of the human reference genome divided by 500 kb chunks (also having the ENCODE blacklist regions excluded); the result is the count matrix showing how many intersections each fragment found in cells in this sample has with each of the 500 kb genomic chunks, sets of the genomic chuncks (index_mapping) and cell barcodes of the sample (barcodes_map)
2. Count the number of fragments in each of the cells in the sample researched

Here is the algorythm of the following analysis:

3. Converting the count matrix with the number of intersections between fragments and chunks to an AnnData object, excluding the cells that have less than 2000 or 4000 fragments from it
4. Filtering doublets via scrublet (functions pp.scrublet and pp.filter_doublets from the Python package SnapATAC2, https://doi.org/10.1038/s41592-023-02139-9 )
5. Excluding 5% of most and less variable chunks
6. Dimensionality reduction using the matrix-free spectral embedding algorithm (function tl.spectral from the Python package SnapATAC2)
7. Converting the data to the 2-dimensional space by UMAP (function tl.umap from the Python package SnapATAC2)
8. Clustering analysis by building a k-nearest neighbour graph and using the Leiden algorithm with resolution 0.5 to identify clusters in the graph (functions pp.knn and tl.leiden from the Python package SnapATAC2)
9. Visualising the data by making the UMAP plot with colloring the plot by Leiden clusterization
10. Annotating cells by number of fragments in them
11. Visualising the data by making the UMAP plot with colloring the plot by number of fragments in cells
10. Making histograms showing distribution of the cells according to the number of fragments in them

Important remark: TSS enrichment scores weren't computed and so filtration of the cells according to this parameter wasn't performed.

The python script should be run this way:

`python python3/scATAC_clusterization_snapATAC_and_figures_making.py ids counts_file input_npz_file barcodes_map index_mapping`

The results are:
1. Umap plots with colloring by the leiden algorythm both for cells containing more than 2000 and 4000 fragments per cell
2. Umap plots with colloring by numbers of fragments in cells both for cells containing more than 2000 and 4000 fragments per cell
3. Histograms of distribution of the cells according to numbers of fragments in cells (for cells containing more than 2000 fragments per cell, from 2000 to 10000 fragments per cell and from 2000 to 10000 fragments per cell but number of cells is in log scale)

