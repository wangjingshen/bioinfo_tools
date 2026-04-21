## Conda env
scFEA_1.1

## Test

script：/SGRNJ03/randd/user/wangjingshen/project/pancreatic_cancer/run_scFEA.py

path: /SGRNJ03/randd/user/wangjingshen/project/pancreatic_cancer/


## Usage

```
python run_scFEA.py --rds data/data_merge.rds --model all --cell_type cell_type --group group  --species human
```
```
help info:
----
--rds                   a rds of seurat object
--group                 name of a variable about group, such as "group"
--cell_type             name of a variable about cell type, such as "cell_type"
--model                 all or split, "all" means using all cells to predict; "split" means using each cell type to predict
--species               human or mouse, set sepecies to get corresponding moduleGene_file, stoichiometry_matrix and cName_file
--n_threads             Optional, number of threads to run scFEA predication, default: 8
--sc_imputation         Optional, whether perform imputation for input (recommend set to True for 10x data), default: False
--cName_file            Optional, the name of compounds
--moduleGene_file       Optional, genes related to metabolic module
--stoichiometry_matrix  Optional, relationship between compounds and modules
```

## Results
#### 01.scFEA_input
  1) input_xxx.csv : scFEA's input data

  2) mapfile.csv : parameters of running scFEA
  
#### 02.scFEA_predication:

1) xxx_module_flux.csv :  The metabolic flux of all metabolic modules in each cell, Specifically, the input (substrate) of the metabolic module decreases corresponding value, and the output (product) of the metabolic module increases the corresponding value;
2) xxx_balance.csv :  The predicted metabolomics;
3) xxx_lossValue.txt :  The loss of the neural network, which converge to a small value during the training;
4) xxx_loss.png :  The loss of the neural network, which converge to a small value during the training;
5) xxx_time.txt :  Run time.

#### 03.scFEA_visualization

1) flux/balance_barplot_xxx.pdf :  Barplot shows difference of flux/balance between group1 and group2. In addtion, the sign of balance represents the accumulate/deplete of the metabolite;

2) flux/balance_xxx.tsv :  Barplot data.
