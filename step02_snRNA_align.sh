#PISA(v0.2;https://github.com/shiquan/PISA/)
#STAR(v2.7.4a)
#sambamba(v.0.7.0)
## Step0 Set parameters
SMAPLE_ID="PG"
DIR1=$PWD
READ1=./fastq1.fq.gz
READ2=./fastq2.fq.gz
INDEX=./Reference/Index/
GTF=./Reference/genomic.gtf

mkdir $DIR1/$SMAPLE_ID ; chmod -R g+s $DIR1/$SMAPLE_ID
mkdir $DIR1/$SMAPLE_ID/temp ; chmod -R g+s $DIR1/$SMAPLE_ID/temp
mkdir $DIR1/$SMAPLE_ID/report ; chmod -R g+s $DIR1/$SMAPLE_ID/report
mkdir $DIR1/$SMAPLE_ID/outs ; chmod -R g+s $DIR1/$SMAPLE_ID/outs

## Step1 PISA_parse
newgrp ST_MCHRI_COHORT
/PISA/bin/PISA parse -t 20 -r 1000000 -f -q 4 -dropN -config /PISA/config/BGI_droplet_scRNA_readStructureV2_T1.json -cbdis $DIR1/$SMAPLE_ID/temp/barcode_counts_raw.txt -report $DIR1/$SMAPLE_ID/report/sequencing_report.csv $READ1 $READ2 -1 $DIR1/$SMAPLE_ID/temp/reads.fq

head -n 100000 $DIR1/$SMAPLE_ID/temp/barcode_counts_raw.txt |cut -f1 > $DIR1/$SMAPLE_ID/temp/barcode_raw_list.txt

## Step2 fastq2bam
/PISA/bin/STAR --outSAMmultNmax 1 --outStd SAM --outSAMunmapped Within --runThreadN 20 --genomeDir $INDEX --readFilesIn $DIR1/$SMAPLE_ID/temp/reads.fq --outFileNamePrefix $DIR1/$SMAPLE_ID/temp/ 1> $DIR1/$SMAPLE_ID/aln.sam && \
/PISA/bin/PISA sam2bam -t 20 -k -o $DIR1/$SMAPLE_ID/temp/aln.bam -report $DIR1/$SMAPLE_ID/report/alignment_report.csv $DIR1/$SMAPLE_ID/aln.sam && rm -f $DIR1/$SMAPLE_ID/aln.sam

## Step3 sortBam
/PISA/bin/sambamba sort -t 20 -o $DIR1/$SMAPLE_ID/temp/sorted.bam $DIR1/$SMAPLE_ID/temp/aln.bam && \
/PISA/bin/PISA anno -gtf $GTF -o $DIR1/$SMAPLE_ID/temp/annotated.bam -report $DIR1/$SMAPLE_ID/report/annotated_report.csv $DIR1/$SMAPLE_ID/temp/sorted.bam && \
/PISA/bin/PISA corr -tag UR -new-tag UB -tags-block CB,GN -@ 20 -cr -o $DIR1/$SMAPLE_ID/outs/final.bam $DIR1/$SMAPLE_ID/temp/annotated.bam

## Step4 cellCount
/PISA/bin/PISA attrcnt -cb CB -tags UB,GN -@ 20 -dedup -o $DIR1/$SMAPLE_ID/temp/cell_stat.txt $DIR1/$SMAPLE_ID/outs/final.bam

## Step5 cellCalling
Rscript /PISA/scripts/scRNA_cell_calling_v1.1.R -i $DIR1/$SMAPLE_ID/temp/cell_stat.txt -o $DIR1/$SMAPLE_ID/report -e 0 -f 0
cp $DIR1/$SMAPLE_ID/report/cell_barcodes.txt $DIR1/$SMAPLE_ID/outs

## Step6 countMatrix
/PISA/bin/PISA count -@ 20 -tag CB -anno-tag GN -umi UB -o $DIR1/$SMAPLE_ID/outs/count_mtx.tsv -list $DIR1/$SMAPLE_ID/outs/cell_barcodes.txt $DIR1/$SMAPLE_ID/outs/final.bam

gzip -f $DIR1/$SMAPLE_ID/outs/count_mtx.tsv
Rscript /PISA/scripts/scRNA_Seurat_clustering.R -i $DIR1/$SMAPLE_ID/outs/count_mtx.tsv.gz -o $DIR1/$SMAPLE_ID/report/

## Step7 report
echo "Name,${SMAPLE_ID}" > /$DIR1/$SMAPLE_ID/report/sample.csv
echo "Original,Null" >>  /$DIR1/$SMAPLE_ID/report/sample.csv
echo "Species,Null" >> /$DIR1/$SMAPLE_ID/report/sample.csv
echo "Sampling time,Null" >> /$DIR1/$SMAPLE_ID/report/sample.csv
echo "Experimental time,Null" >> /$DIR1/$SMAPLE_ID/report/sample.csv
/JAVA/jdk-12.0.1/bin/java -jar /PISA/scripts/idrop-0.0.1.jar $DIR1/$SMAPLE_ID/report $SMAPLE_ID

