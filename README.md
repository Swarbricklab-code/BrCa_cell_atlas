# BrCa_cell_atlas
This repository contains code related to data processing and downstream analysis associated with the study "[A single-cell and spatially resolved atlas of human breast cancers](https://www.nature.com/articles/s41588-021-00911-1)" at Nature Genetics. 

# Data Availability
## Processed & filtered scRNA-Seq data and CITE data
All processed scRNA-seq data and is available for in-browser exploration and download through the [Broad Single-Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1039). 
The same data is also available via [CELLxGENE](https://cellxgene.cziscience.com/collections/65db5560-7aeb-4c66-b150-5bd914480eb8).
CITE data can also be downloaded from this portal, although it cannot be explored interactively at this point.
Processed data is also available on GEO ([GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)).

## Raw scRNA-Seq data
Raw bulk RNA-seq data and scRNA-Seq data from this study has been deposited in the European Genome-Phenome Archive (EGA), which is hosted by the EBI and the CRG, under the accession code [EGAS00001005173](https://ega-archive.org/studies/EGAS00001005173). 

## Spatial transcriptomics data
All ST data from this study is now available on [CELLxGENE](https://cellxgene.cziscience.com/collections/65db5560-7aeb-4c66-b150-5bd914480eb8).
Intermediate files are available from the Zenodo data repository (DOI: 10.5281/zenodo.4739739 - https://zenodo.org/record/4739739). ST data from the Andersson et al. study can be downloaded from the Zenodo data repository (DOI: 10.5281/zenodo.3957257 - https://zenodo.org/record/4751624#.Yw06i3ZByMI). For deconvolution of Visium data, please see the [stereoscope method](https://github.com/almaan/stereoscope) from the [Andersson et al (2020) study](https://www.nature.com/articles/s42003-020-01247-y).
See [this issue thread](https://github.com/Swarbricklab-code/BrCa_cell_atlas/issues/3) for a discussion on how to load the spatial data from Zenodo into Seurat.

# Contacts
All other relevant data and analysis are available from the authors upon request. For further enquires, please either raise an issue via GitHub or email John Reeves (Data Manager - j.reeves(at)garvan.org.au), Sunny Wu (Lead Author - s.wu(at)garvan.org.au) or Alexander Swarbrick (Lab Head - a.swarbrick(at)garvan.org.au).
