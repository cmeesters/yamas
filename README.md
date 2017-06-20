[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# yamas

YAMAS is a software for meta-analysis of genome wide association studies.

A main purpose of imputation prior to metaanalysis (MA) is to unify the available marker panels of genome wide analyses (GWAS) and to avoid loss of SNPs that are present in one study but not in another. With YAMAS we present a MA software that avoids such loss without the need to impute data. By using reference data from the HAPMAP and 1000 genomes projects users are enabled to analyse all SNPs that are present in at least one of the experimental marker panels: LD-information is used to find substitute makers for those missing.

## Notes to this repository
* This repository is a clone of [the original source](http://yamas.meb.uni-bonn.de/), which was created using svn.
* Please see the reference below and acknowledge the original authors and cite the reference, if using this software.

## Current features include:

    meta-analysis marker by marker with fixed and random effects
    meta-analysis per LD-block
    filling up untyped or missing markers with proxy markers prior to meta-analysis – a reference file containing r2-values from HAPMAP/1k is provided
    parallel reading of files
    YAMAS is able to parse any tabulate file format
    some algorithms are parallelized


## Citation:
Meesters C, Leber M, Herold C, Angisch M, Mattheisen M, Drichel D, Lacour A, Becker T (2012). Quick, "imputation-free" meta-analysis with proxy-SNPs. BMC Bioinformatics. 2012 Sep 12;13:231. doi: 10.1186/1471-2105-13-231.
