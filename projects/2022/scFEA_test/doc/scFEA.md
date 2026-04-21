Here, we use scFEA to estimate cell-wise metabolic flux using single cell RNA-seq data.

### Conda env
scFEA_1.1  (/SGRNJ/Public/Software/conda_env/scFEA_1.1)

### Path
  scFEA :  /SGRNJ03/randd/user/wangjingshen/bio_soft/scFEA_v1.1.2/

  test : /SGRNJ03/randd/user/wangjingshen/project/metabolism/

### Input 
scRNA-seq expression matrix (raw count or normalised counts)

### Run
1.run with default parameters
```
./run_scFEA.sh -t ../data/Melissa_full.csv
```

2.run with mouse data (The default model is human)
```
./run_scFEA.sh -m mouse -t ../data/mouse_example_data.csv -i True
```

3.run with interested metabolic processes(such as the "glutaminolysis" provided by the author)
```
./run_scFEA.sh -t ../data/Melissa_full.csv \
  -f ../../../bio_soft/scFEA_v1.1.2/data/module_gene_glutaminolysis1_m23.csv \
  -s ../../../bio_soft/scFEA_v1.1.2/data/cmMat_glutaminolysis1_c17_m23.csv \
  -c ../../../bio_soft/scFEA_v1.1.2/data/cName_glutaminolysis1_c17_m23.csv \
  -o ../output/output_test/ 
```

4.get help information 
```
$ ./run_scFEA.sh -h
```
```
usage: run_scFEA.sh [--test_file ] [--modle ] [--output_path ] [--sc_imputation ] [--n_threads ] [--help]
 -t, --test_file
         A single cell expression matrix, and it can be raw counts or normalised counts
 -m, --model
         Optional, predicting metabolism to human or mouse. default: human
 -f, --moduleGene_file
         Optional, the table contains genes for each module. default: module_gene_m168.csv(human)
 -s, --stoichiometry_matrix
         Optional, the table describes relationship between compounds and modules. default: cmMat_c70_m168.csv(human)
 -c, --cName_file
         Optional, the name of compounds. default: cName_c70_m168.csv(human)
 -o, --output_path
         Optional, the data directory for result. default: ../output/
 -i, --imputation
         Optional, whether perform imputation for SC dataset (recommend set to <True> for 10x data). default: False
 -p, --n_threads
         Optional, number of threads used for artificial neural network calculations. default: 8
 -h, --help
         Optional, show help message and exit.
```

### Output

Here, I take the result of Melissa in output as an example:

1)Melissa_full_module_flux.csv

&nbsp;  The metabolic flux of all metabolic modules in each cell, Specifically, the input (substrate) of the metabolic module decreases corresponding value, and the output (product) of the metabolic module increases the corresponding value;

2)Melissa_full_balance.csv

&nbsp;  The predicted amount of all metabolites in each cell (metabolic imbalance), this value is calculated by subtracting inputs (substrates) of all metabolic module from outputs (products) of all metabolic module;

3)Melissa_full_lossValue.txt and &nbsp;4)Melissa_full_loss.png

&nbsp;  The loss of the neural network, which converge to a small value during the training of the scFEA model;

5)Melissa_full_time.txt

&nbsp; Run time.

