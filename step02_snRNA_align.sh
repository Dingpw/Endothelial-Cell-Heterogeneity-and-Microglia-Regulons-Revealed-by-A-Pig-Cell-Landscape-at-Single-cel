## You should download and install these three tools before running this step:
## Cromwell release 35 (https://github.com/broadinstitute/cromwell/)
## DNBelab_C_Series_HT_scRNA-analysis-software (https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)
## PISA version 0.5 (https://github.com/shiquan/PISA/)
## R version 3.5.2 (https://www.r-project.org/). R packages: “ggplot2”, “getopt”, “data.table”, “cowplot”, “DropletUtils”, and “Seurat”.

## "Droplet_scRNA_pipeline_html.wdl" file was stored in pipeline "DNBelab_C_Series_HT_scRNA-analysis-software"
java -jar tools/DNBelab_C_Series_HT_scRNA-analysis-software/cromwell-35.jar run -i config.json pipelines/Droplet_scRNA_pipeline_html.wdl
