#! /usr/bin/env bash

echo "This script is used for transfer scanpy h5ad to Seurat rds."
read -p "Enter the path of h5ad: " h5ad_path
read -p "Whether to tranfer the Seurat V4 object to Seurat V3 object (yes or no):" v4_to_v3
read -p "Prefix of output files:" prefix
read -p "Path storing ouput files:" outdir
read -p "Wether perform normalize(yes or no):" norm 

outdir=`echo ${outdir}|sed 's/\"//g'`
if [ ! -e ${outdir} ];then
	mkdir ${outdir}
else
	echo "outdir existed"
fi

h5ad_path=`echo ${h5ad_path}|sed 's/\"//g'`

echo "h5ad_path: ${h5ad_path}
wether transfer v4 to v3: ${v4_to_v3}
prefix: ${prefix}
outdir: ${outdir}" > ${outdir}/logfile

cat ${outdir}/logfile


source activate sceasy
echo "#! /usr/bin/env Rscript
library(sceasy)
convertFormat('${h5ad_path}',from = 'anndata',to = 'seurat',main_layer = 'counts',outFile = '${outdir}/${prefix}_h5ad_to_rds.rds')
data <- readRDS('${outdir}/${prefix}_h5ad_to_rds.rds')
print(colnames(data@meta.data))
data@assays$RNA@key <- "rna_"   # fix
saveRDS(data,'${outdir}/${prefix}_h5ad_to_rds.rds')
print('convert finished')
" > ./convert.R
echo "R --no-save < convert.R" > ./run.sh
if [ ${v4_to_v3} == yes ];then
#	source activate Seuratv3
#	mkdir ${outdir}/script
#	cd ${outdir}/script
	echo "/SGRNJ06/randd/USER/wangjingshen/script/data_trans/script/tran_v4_to_v3.R \
	--v4rds ${outdir}/${prefix}_h5ad_to_rds.rds \
	--prefix ${prefix} \
	--cluster_col cluster \
	--scale_data no \
	--normalize ${norm} \
	--outdir ${outdir}" > ./v4_to_v3.sh
fi

echo `[ ${v4_to_v3} == no ]`
if [ ${v4_to_v3} == no ];then
	echo "source activate sceasy 
bash run.sh" > trans_h5ad_to_seuratv4.sh
else
	echo "source activate sceasy
bash run.sh
source activate Seuratv3
bash v4_to_v3.sh" > trans_h5ad_to_seuratv3.sh
fi
