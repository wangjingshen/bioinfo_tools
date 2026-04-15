Here, we use calculate_gene_fraction_reads_in_cells.py to calculate the fraction reads in cells for each gene.

## Conda env
scFEA_1.1

## Test
script：/SGRNJ06/randd/USER/wangjingshen/script/2022/gene_fraction_reads_in_cells/script/calculate_gene_fraction_reads_in_cells.py

path:  /SGRNJ06/randd/USER/wangjingshen/script/2022/gene_fraction_reads_in_cells/test/

## Usage
```
python calculate_gene_fraction_reads_in_cells.py \
  --counts /SGRNJ03/randd/user/zhouyiqi/project/full_len_vdj/P-S1TDE/05.count/P-S1TDE_counts.txt \
  --count_detail /SGRNJ03/randd/user/zhouyiqi/project/full_len_vdj/P-S1TDE/05.count/P-S1TDE_count_detail.txt \
  --outdir ../result/
```
```
help info:
----
  -h, --help                     show this help message and exit
  --counts COUNTS                05.count/xxx_counts.txt, a output of CeleScope
  --count_detail COUNT_DETAIL    05.count/xxx_count_detail.txt, a output of CeleScope
  --outdir OUTDIR                output dir
```

## Result
gene_fraction_reads_in_cells.tsv (fraction_reads_in_cells were kept to 4 decimal places)
```
geneID  fraction_reads_in_cells
ENSG00000000003 0.6792
ENSG00000000419 0.7058
ENSG00000000457 0.7201
ENSG00000000460 0.6803
ENSG00000000938 0.5225
...
```
