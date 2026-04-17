import pandas as pd
import argparse
import subprocess



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
    return(cmd)


    #tab = '\t'
    #cmd = f"""awk -F "\t" '$3=="gene" {
    #    match($9, /gene_id "([^"]+)"/, gid);
    #    match($9, /gene_name "([^"]+)"/, gname);
    #    match($9, /gene_biotype "([^"]+)"/, gbio);
    #    print gid[1], gname[1], gbio[1]
    #}' {gtf} | sed 's/ /{tab}/g' | sort > {sample}_biotype.tsv """



