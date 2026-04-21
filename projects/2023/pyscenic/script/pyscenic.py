import os
import argparse
import subprocess
from utils import h5ad2loom
from utils import mtx2loom


class Pyscenic:
    def __init__(self, args):
        self.input = args.input
        self.mode = args.mode
        self.species = args.species
        self.threads = args.threads
        self.subset = args.subset 
        self.outdir = args.outdir
        if(self.species=='human'):
            self.tf_files = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/tf/allTFs_hg38.txt'
            self.rank_feature = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/cistarget/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/cistarget/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
            self.tbl = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/motif/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
        if(self.species=='mouse'):
            self.tf_files = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/tf/allTFs_mm.txt'
            self.rank_feature = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/cistarget/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/cistarget/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
            self.tbl = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'
        if(self.species=='fly'):
            self.tf_files = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/tf/allTFs_dmel.txt'
            self.rank_feature = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/cistarget/dm6_v10_clust.genes_vs_motifs.rankings.feather'
            self.tbl = '/SGRNJ06/randd/USER/wangjingshen/script/pyscenic/data/motif/motifs-v10nr_clust-nr.flybase-m0.001-o0.0.tbl'

        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}')
    
    def get_loom(self):
        '''
        1) get loom from annotated h5ad (recommended); 
        2) loom;
        3) mtx to loom (Not recommended. It is recommended to use mtx_to_h5ad_loom to generate h5ad and loom);
        '''
        if(self.mode == 'h5ad'):
            h5ad2loom(self.input, outdir = self.outdir)
            self.loom = f'{self.outdir}/pyscenic_input.loom'
        if(self.mode == 'loom'):
            self.loom = self.input
        if(self.mode == 'mtx'):
            mtx2loom(self.input, outdir = self.outdir)
            self.loom = f'{self.outdir}/pyscenic_input.loom'

    def run_pyscenic(self):
        cmd1 = f"/SGRNJ/Public/Software/conda_env/pyscenic_env/bin/pyscenic grn {self.loom} {self.tf_files} -o {self.outdir}/pyscenic_adjacencies.csv --num_workers {self.threads}"
        cmd2 = f"/SGRNJ/Public/Software/conda_env/pyscenic_env/bin/pyscenic ctx {self.outdir}/pyscenic_adjacencies.csv {self.rank_feature} --annotations_fname {self.tbl} --expression_mtx_fname {self.loom} --output {self.outdir}/pyscenic_regulons.csv --mask_dropouts --num_workers {self.threads}"
        cmd3 = f"/SGRNJ/Public/Software/conda_env/pyscenic_env/bin/pyscenic aucell {self.loom} {self.outdir}/pyscenic_regulons.csv --output {self.outdir}/pyscenic_output.loom --num_workers {self.threads}"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)

    def reg2cytoscape(self):
        if(self.subset == None):
            cmd1 = f"/SGRNJ/Public/Software/conda_env/jsr4.1/bin/Rscript /SGRNJ06/randd/USER/wangjingshen/script/pyscenic/script/regulons_cytoscape.R --reg {self.outdir}/pyscenic_regulons.csv --outdir {self.outdir}"
        else:
            cmd1 = f"/SGRNJ/Public/Software/conda_env/jsr4.1/bin/Rscript /SGRNJ06/randd/USER/wangjingshen/script/pyscenic/script/regulons_cytoscape.R --reg {self.outdir}/pyscenic_regulons.csv --subset {self.subset} --outdir {self.outdir}"
        subprocess.check_call(cmd1, shell = True)

    def run(self):
        self.get_loom()
        self.run_pyscenic()
        self.reg2cytoscape()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--input', help='input file', required=True)
    parsers.add_argument('--mode', help='mode, h5ad(default) or loom or mtx')
    parsers.add_argument('--species', help='species, human or mouse or fly', required=True)
    parsers.add_argument('--threads', help='threads', required=True)
    parsers.add_argument('--subset', help='subset')
    parsers.add_argument('--outdir', help='outdir', required=True)
    args = parsers.parse_args()
    if(args.mode == None):
        args.mode = 'h5ad'
    runner = Pyscenic(args) 
    runner.run()

if __name__ == '__main__':
    main()





