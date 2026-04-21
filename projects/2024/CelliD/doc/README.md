Here, we use CelliD to automatically annotate cell type.

## Conda env
jsr4.1

## Usage
run_CelliD.R
```
--rds          seurat rds file
--resolution   Optional, resolution, Default: 'seurat_clusters' in rds
--species      Optional, species, Hs(Default, human) or Mm(mouse)
--ref          Optional, markers reference, raw(Default, PanglaoDB_markers) or update(markers from data scienece)
--mode         analysis mode, all(markers of all organ) or single(markers of single organ) or common(markers of single organ and common celltype)
--organ        name of organ, raw reference: Adrenal glands, Blood, Bone, Brain, Connective tissue, Embryo, Epithelium, Eye, GI tract, Heart, Immune system, Kidney,
                                             Liver, Lungs, Mammary gland, Olfactory system, Oral cavity, Pancreas, Placenta, Reproductive, Skeletal muscle, Skin, 
                                             Smooth muscle, Thymus, Thyroid, Urinary bladder, Vasculature, Zygote
                              update reference: BM_PBMC
--nfeatures    Optional, integer of top n features to consider for hypergeometric test, Default: 200, larger number(for example, 700) may be suitable for update ref
--outdir       Optional, outdir, Default: ./CelliD_outdir
```

## Test path

1)run single samples:

/SGRNJ06/randd/USER/wangjingshen/project/CelliD/test/test_single_sample/pancreatic_cancer_mode_common/run.sh

2)run multiple samples:

/SGRNJ06/randd/USER/wangjingshen/project/CelliD/test/test_multi_samples/run.sh

## Input

an example of mapfile for multiple samples 

```
rds     resolution      species mode    organ   ref     nfeatures       outdir
/SGRNJ06/randd/USER/wangjingshen/project/pancreatic_cancer/data/data_merge_new.rds      0.6     Hs      single  Pancreas        raw     200     mode_single
/SGRNJ06/randd/USER/wangjingshen/project/pancreatic_cancer/data/data_merge_new.rds      0.6     Hs      all     Pancreas        raw     200     mode_all
/SGRNJ06/randd/USER/wangjingshen/project/pancreatic_cancer/data/data_merge_new.rds      0.6     Hs      common  Pancreas        raw     200     mode_common
```

## Output

1)CelliD_prediction.pdf: 

figure of cell annotation 
  
2)CelliD_prediction.txt: 

cell annotation 
  
3)CelliD_prediction_freq_each_cluster.txt: 

prediction frequency in each cluster
