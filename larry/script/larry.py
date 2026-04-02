import gzip
import time 
import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

larry_prefix = 'GTTGCTAGGAGAGACCATATG'
N_READS = 10
N_UMIS = 3
N_HAMMING = 3

def is_valid(bc):
    return (bc[4:6]=='TG'
            and bc[10:12]=='CA' 
            and bc[16:18]=='AC' 
            and bc[22:24]=='GA' 
            and bc[28:30]=='GT'
            and bc[34:36]=='AG')

def in_filtered_list(cell_bc, filtered_cell_barcodes):
    num_Ns = sum([c=='N' for c in cell_bc])
    if num_Ns > 1: return False
    elif num_Ns == 1: return np.any([cell_bc.replace('N',c) in filtered_cell_barcodes for c in 'ACTG'])
    else: return cell_bc in filtered_cell_barcodes


def larry(spname, library_id, outdir):
    sample_paths = {
        f'{spname}': {
            'filtered_cell_barcodes_path': f'{spname}/input/matrix_10X/barcodes.tsv.gz',
            'larry_barcode_fastq_paths':[f'{spname}/input/{library_id}'"_{}.fastq.gz"]
        }
    }


    counts = {}

    for sample,paths in sample_paths.items():
        filtered_cell_barcodes = gzip.open(paths['filtered_cell_barcodes_path']).read().decode('utf-8').split('\n')
        filtered_cell_barcodes_set = set(filtered_cell_barcodes)

        for fastq_path in paths['larry_barcode_fastq_paths']:
            R1 = gzip.open(fastq_path.format('R1'))
            R2 = gzip.open(fastq_path.format('R2'))
            counter = 0
            start_time = time.time()
            while True:
                counter += 1
                if counter % 1000000 == 0: print(fastq_path+ ': Processed {} lines in {} seconds'.format(counter, time.time()-start_time))
                try:
                    r1_line = R1.readline().decode('utf-8')
                    r2_line = R2.readline().decode('utf-8')
                except:
                    print('ERROR extracting {}'.format(fastq_path))
                    break
                if r2_line == '': break
                if r2_line[0] in '@+': continue
                if larry_prefix in r2_line:
                    larry_bc = r2_line.split(larry_prefix)[1][:40]
                    cell_bc =  r1_line[0:9] + "_" + r1_line[15:24] + "_" + r1_line[30:39]         # GEXSCOPE-V3 pattern C9L6C9L6C9L1U12   # 10X r1_line[:16]+'-1'
                    umi = r1_line[40:52]         # 10X r1_line[16:24]
                    if is_valid(larry_bc) and in_filtered_list(cell_bc, filtered_cell_barcodes):
                        combo = (sample, cell_bc, umi, larry_bc)
                        if combo in counts:
                            counts[combo] += 1
                        else:
                            counts[combo] = 1

    # Set parameters for clonal analsysis
    num_reads = [v for k,v in counts.items()]
    plt.hist(np.log(num_reads)/np.log(10), bins=50)
    plt.axvline(np.log(N_READS)/np.log(10),c='k')
    plt.xticks(range(5),np.logspace(0,4,5))
    plt.text(np.log(N_READS)/np.log(10)*1.1,10**6,'N_READS cutoff', fontsize=12)
    plt.yscale('log')
    plt.savefig(f'{outdir}/N_READS_cutoff.png')
    plt.close()

    counts_filtered = {k:v for k,v in counts.items() if v >= N_READS}
    print('Retaining '+repr(len(counts_filtered))+ ' out of '+repr(len(counts))+' (Sample,Cell-BC,UMI,GFP-BC) combinations')


    # Collapse GFP-BCs by hamming distance
    def hamming(bc1,bc2): return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])

    all_gfp_bcs = sorted(set([k[3] for k in counts_filtered]))
    good_gfp_bcs = []
    bc_map = {}
    for i,bc1 in enumerate(all_gfp_bcs):
        if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(all_gfp_bcs))+' barcodes')
        mapped = False
        for bc2 in good_gfp_bcs:
            if hamming(bc1,bc2) <= N_HAMMING:
                mapped = True
                bc_map[bc1] = bc2
                break
        if not mapped:
            good_gfp_bcs.append(bc1)

    print('\nCollapsed '+repr(len(bc_map))+' barcodes')
    for bc in good_gfp_bcs: bc_map[bc] = bc


    # Filter GFP-barcodes by UMI
    cell_data = {}
    for sample,paths in sample_paths.items():
        filtered_cell_barcodes = gzip.open(paths['filtered_cell_barcodes_path']).read().decode('utf-8').split('\n')
        for cell_bc in filtered_cell_barcodes:
            cell_data[(sample,cell_bc)] = {}

    for sample,cell_bc,umi,larry_bc in counts_filtered.keys():
        if (sample,cell_bc) in cell_data:
            if not larry_bc in cell_data[(sample,cell_bc)]:
                cell_data[(sample,cell_bc)][larry_bc] = 0
            cell_data[(sample,cell_bc)][larry_bc] += 1

    num_cells_with_barcode = np.zeros(20)
    for larry_bc_counts in cell_data.values():
        if len(larry_bc_counts)>0:
            num_cells_with_barcode[:np.min([20,np.max(list(larry_bc_counts.values()))])] += 1
    efficiency = num_cells_with_barcode / len(cell_data)
    plt.plot(range(20),efficiency)
    plt.plot([N_UMIS,N_UMIS],[np.min(efficiency),np.max(efficiency)],'-k',linewidth=2)
    plt.text(N_UMIS*1.1,np.max(efficiency)*.95,'UMI cutoff',fontsize=14)
    plt.savefig(f'{outdir}/UMI_cutoff.png')
    plt.close()

    final_BCs = {}
    for k,larry_bc_counts in cell_data.items():
        final_BCs[k] = '-'.join(sorted([k for k,v in larry_bc_counts.items() if v >= N_UMIS]))
    print('\nFinal annotation has '+repr(len(set(final_BCs.values())))+' clones in '+repr(len([k for k,v in final_BCs.items() if len(v)>0]))+' cells')


    output = []
    for sample,paths in sample_paths.items():
        filtered_cell_barcodes = gzip.open(paths['filtered_cell_barcodes_path']).read().decode('utf-8').split('\n')
        for cell_bc in filtered_cell_barcodes:
            output.append(sample+','+cell_bc+','+final_BCs[(sample,cell_bc)])
    open(f'{outdir}/{spname}.larry_clones.csv','w').write('\n'.join(output))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--spname', help='spname', required = True)
    parser.add_argument('--library_id', help='library id', required = True)
    parser.add_argument('--outdir', help='outdir', required = True)
    args = parser.parse_args()

    runner = lappry(args.spname, args.library_id, args.outdir)
    runner.run()


if __name__ == '__main__':
    main()