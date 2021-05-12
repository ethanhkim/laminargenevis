## LaminaRGeneVis: a tool to visualize a tool to visualize gene expression across the laminar architecture of the human neocortex
[<img src="https://img.shields.io/badge/Supported%20by-%20CONP%2FPCNO-red?link=https://conp.ca/">](https://conp.ca/)
<!---[![DOI:<your number>](http://img.shields.io/badge/DOI-<your number>-<colour hexcode>.svg)](<doi link>)--->

Welcome to LaminaRGeneVis! 
This [web application](https://ethanhkim.shinyapps.io/laminargenevis/) allows you to examine layer-specific gene expression across the cortex and determine layer annotations.

## App Workflow:
Here's a quick overview of the app! Multiple datasets of RNA-seq expression (described below) have been standardized for ease of comparison as shown in (A). Users have the option of choosing to examine either Single Gene or Multiple Genes in the HGNC gene symbol format. Choosing Single Gene will give you a plot similar to (B), and choosing Multiple Genes will give you multiple plots, including the plot shown in (C).

![Workflow Figure](https://github.com/ethanhkim/laminargenevis/blob/master/www/pageFigure.png)

## Source Data:
The data for this application has been sourced from these following studies and institutions:

1. [He et al. (2017): ](https://pubmed.ncbi.nlm.nih.gov/28414332/) 
  * This study assayed the whole genome using high-throughput RNA-seq in samples from the DLPFC.
2. [Maynard et al. (2021): ](https://www.nature.com/articles/s41593-020-00787-0) 
  * This study assayed the whole genome through the Visium Platform (10X Genomics) in samples from the DLPFC.
3. [Allen Institute for Brain Science (AIBS): Cell-Type Database](https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq)
  * This dataset contains multiple cortical regions and assays the whole genome across roughly 49,000 single-cell nuclei.


## Data processing and availability:
* Data processing scripts are available under [/raw_data_processing](https://github.com/ethanhkim/laminargenevis/tree/master/raw_data_processing). 
* For processing the raw files for He et al, the scripts are available [here](https://github.com/derekhoward/he_seq). 
* Processed data are available under [/data/processed](https://github.com/ethanhkim/laminargenevis/tree/master/data/processed). 

## Support
This work was supported in part by funding provided NSERC, and by [Brain Canada](https://braincanada.ca/), in partnership with [Health Canada](https://www.canada.ca/en/health-canada.html), for the [Canadian Open Neuroscience Platform initiative](https://conp.ca/).


[<img src=”https://conp.ca/wp-content/uploads/elementor/thumbs/logo-2-o5e91uhlc138896v1b03o2dg8nwvxyv3pssdrkjv5a.png”>](https://conp.ca/)