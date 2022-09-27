## CellRanger version 3.0.2 (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

/Software/cellranger-3.0.2_R/cellranger count \
                  --id=Pig-Adipose-S \
                  --transcriptome=/dellfsqd2/ST_LBI/USER/liangxue/reference/refdata-cellranger-pig/susScr11_genome \
                  --fastqs=/dellfsqd2/ST_LBI/USER/liangxue/project/PigAtlas/Pig-Adipose-S/input/ \
                  --sample=Adipose-S \
                  --expect-cells=10000 \ 
                  --jobmode=local \
                  --localcores 10 \
                  --localmem 50
