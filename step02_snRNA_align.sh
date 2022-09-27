## You should download and install these three tools before running this step:
## Cromwell release 35 (https://github.com/broadinstitute/cromwell/releases/tag/35)
## DNBelab_C_Series_HT_scRNA-analysis-software (https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software)
## PISA v0.5 (https://github.com/shiquan/PISA/)

## "Droplet_scRNA_pipeline_html.wdl" file was stored in pipeline "DNBelab_C_Series_HT_scRNA-analysis-software"
java -jar tools/DNBelab_C_Series_HT_scRNA-analysis-software/cromwell-35.jar run -i config.json pipelines/Droplet_scRNA_pipeline_html.wdl
