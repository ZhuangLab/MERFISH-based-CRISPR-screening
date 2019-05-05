# MERFISH-based_CRISPR_screening_in_mammalian_cells

This project contains code for image and sequence analysis for pooled CRISPR screening using MERFISH for genotype identification in individual cells as described in Imaging-based pooled CRISPR screening reveals regulators of lncRNA localization.

Please check github.com/ZhuangLab for the latest version of this code.
All function libraries can be found in github.com/ZhuangLab

-----Matlab image analysis-----

This analysis software is used for identifying segmenting images into cells, extracting barcodes from each cell, and analyzing cell phenotypes. This analysis software is written in matlabs, where the path is 

------Matlab sequence analysis------

The Matlab sequence analysis for building the barcode to genetic variant lookup table is comprised of three parts in this scripts. 1. parse reads from miSeq to identify the barcode and Protospacers within the reads 2. find unique protospacer-barcode reads within the library 3. Further identify the unique reads considering sequencing errors within the unique protospacer-barcode reads to generate the final lookup table.


Code authors

Tian Lu [lutian125@gmail.com]
