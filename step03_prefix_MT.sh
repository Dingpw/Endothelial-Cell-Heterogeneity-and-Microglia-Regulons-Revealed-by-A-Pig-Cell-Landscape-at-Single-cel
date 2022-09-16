# Add“MT-” prefix for the pig mitochondrial genes in the “features.csv”(scRNA-seq) and "matrix.mtx"(snRNA-seq) file of each sample manually before Seurat.
## scRNA-seq
sed -i 's/\<ND2\>/MT-ND2/g' features_MT.tsv 
sed -i 's/\<COX1\>/MT-COX1/g' features_MT.tsv 
sed -i 's/\<COX2\>/MT-COX2/g' features_MT.tsv 
sed -i 's/\<ATP8\>/MT-ATP8/g' features_MT.tsv 
sed -i 's/\<ATP6\>/MT-ATP6/g' features_MT.tsv 
sed -i 's/\<COX3\>/MT-COX3/g' features_MT.tsv 
sed -i 's/\<ND3\>/MT-ND3/g' features_MT.tsv 
sed -i 's/\<ND4L\>/MT-ND4L/g' features_MT.tsv 
sed -i 's/\<ND4\>/MT-ND4/g' features_MT.tsv 
sed -i 's/\<ND5\>/MT-ND5/g' features_MT.tsv 
sed -i 's/\<ND6\>/MT-ND6/g' features_MT.tsv 
sed -i 's/\<CYTB\>/MT-CYTB/g' features_MT.tsv

## snRNA-seq
sed -i 's/\<ND2\>/MT-ND2/g' matrix_MT.mtx 
sed -i 's/\<COX1\>/MT-COX1/g' matrix_MT.mtx 
sed -i 's/\<COX2\>/MT-COX2/g' matrix_MT.mtx  
sed -i 's/\<ATP8\>/MT-ATP8/g' matrix_MT.mtx  
sed -i 's/\<ATP6\>/MT-ATP6/g' matrix_MT.mtx  
sed -i 's/\<COX3\>/MT-COX3/g' matrix_MT.mtx  
sed -i 's/\<ND3\>/MT-ND3/g' matrix_MT.mtx  
sed -i 's/\<ND4L\>/MT-ND4L/g' matrix_MT.mtx  
sed -i 's/\<ND4\>/MT-ND4/g' matrix_MT.mtx  
sed -i 's/\<ND5\>/MT-ND5/g' matrix_MT.mtx  
sed -i 's/\<ND6\>/MT-ND6/g' matrix_MT.mtx  
sed -i 's/\<CYTB\>/MT-CYTB/g' matrix_MT.mtx 
