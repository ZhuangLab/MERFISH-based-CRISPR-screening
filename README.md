# MERFISH-based_CRISPR_screening_in_mammalian_cells

This project contains code for image and sequence analysis for pooled CRISPR screening using MERFISH for genotype identification in individual cells as described in Imaging-based pooled CRISPR screening reveals regulators of lncRNA localization.

Please check github.com/ZhuangLab for the latest version of this code.

Example data can be found at zhuang.harvard.edu/merfish.html.

-----Matlab image analysis-----

This analysis software is used for segmenting images into cells, extracting barcodes from each cell, and analyzing cell phenotypes. This analysis software is written in an object-oriented manner and found in the core directory. Analysis is controlled through core\experiment.py, which takes the following parameters:

------Matlab sequence analysis------

The Matlab sequence analysis for building the barcode to genetic variant lookup table is comprised of three scripts


Code authors

Tian Lu [lutian125@gmail.com]
