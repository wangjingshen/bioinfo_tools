# import dependencies
import os
import argparse
import scanpy as sc


class H5ad_anno:
    def __init__(self, args):
        self.h5ad = args.h5ad
        self.anno = args.anno
        self.mode = args.mode
        self.outdir = args.outdir
        self.name = args.name
    
    def annotation(self):
        adata = sc.read_h5ad(self.h5ad)
        if self.mode=='cmd':
            dict_1 ={}
            for i,anno_i in zip(adata.obs['cluster'].cat.categories, self.anno.split(',')):
                dict_1[i] = anno_i
        
        if self.mode == 'csv':
            with open(self.anno) as fin:
                anno = [line.rstrip() for line in fin]
                dict_1 = dict([x.split(',') for x in anno])
        
        adata.obs['cell_type'] = adata.obs['cluster'].map(dict_1)
        adata.write(f'{self.outdir}/{self.name}.h5ad')
    

    def run(self):
        self.annotation()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--h5ad', help = 'h5ad', required = True)
    parsers.add_argument('--anno', help = 'anno',required = True)
    parsers.add_argument('--mode', help = 'mode',required = True)
    parsers.add_argument('--outdir', help = 'outdir', required = True)
    parsers.add_argument('--name', help = 'name', required = True)
    args = parsers.parse_args()

    runner = H5ad_anno(args) 
    runner.run()

if __name__ == '__main__':
    main()
