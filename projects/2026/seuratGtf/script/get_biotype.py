import subprocess
import pandas as pd
import argparse
import os
import sys
from pathlib import Path
import time

def insert_sys_path(levels_up):
    '''
    insert bioinfo_tools to sys.path
    '''
    root = Path(__file__).resolve()
    for _ in range(levels_up):
        root = root.parent
    if not (root / "utils").exists():
        raise FileNotFoundError(f"utils not found in {root}.")
    sys.path.insert(0, str(root))
    return(root)

bioinfo_root = insert_sys_path(5)  # Top 5 parent directories of current script (bioinfo_tools)
pipeline_root = Path(__file__).resolve().parents[1]

from utils.utils import execute_cmd, timer


def get_biotype_from_gtf(gtf, outdir):
    '''
    some ens_id no gene_name, use ens_id
    # grep -A 2 -B 2 ENSMUSG00002076890 mus_110.tsv
    '''
    tab = '\t'
    cmd = """awk -F "\\t" '$3=="gene" {{
        match($9, /gene_id "([^"]+)"/, gid);
        if(match($9, /gene_name "([^"]+)"/, gname)){{
            gene_name = gname[1];
        }}else{{
            gene_name = gid[1];
        }}
        match($9, /gene_biotype "([^"]+)"/, gbio);
        print gid[1], gene_name, gbio[1]
        }}' {gtf} | sed 's/ /{tab}/g' | sort > {outdir}/gtf_biotype.tsv
    """.format(gtf=gtf, tab=tab, outdir=outdir)
    execute_cmd(cmd)



    #tab = '\t'
    #cmd = f"""awk -F "\t" '$3=="gene" {
    #    match($9, /gene_id "([^"]+)"/, gid);
    #    match($9, /gene_name "([^"]+)"/, gname);
    #    match($9, /gene_biotype "([^"]+)"/, gbio);
    #    print gid[1], gname[1], gbio[1]
    #}' {gtf} | sed 's/ /{tab}/g' | sort > {sample}_biotype.tsv """

def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--gtf', help='gtf', required=True)
    parsers.add_argument('--outdir', default = "outdir", help='outdir, default: outdir')

    args = parsers.parse_args()
    get_biotype_from_gtf(args.gtf, args.outdir) 

if __name__ == '__main__':
    main()


