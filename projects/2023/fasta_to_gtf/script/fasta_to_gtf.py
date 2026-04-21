import os
import argparse
import subprocess

class Fasta_to_gtf:
    def __init__(self, args):
        self.genome = args.genome
        self.fasta = args.fasta
        self.prefix = args.prefix
        self.star_version = args.star_version
        self.update = args.update

        if not os.path.exists("star"):
            os.system(f"mkdir -p star")

        if self.star_version == "new":
            self.star = "/SGRNJ/Public/Software/conda_env/celescope2.1.0/bin/STAR-avx2"
        else:
            self.star = "/Public/Software/miniconda2/envs/old/bin/STAR"


    def run(self):
        '''
        cmd1 change dir, so fasta need absolute Path
        '''
        cmd1 = f"{self.star} --runThreadN 4 \
                    --genomeDir {self.genome} \
                    --readFilesIn {self.fasta} \
                    --outSAMtype BAM SortedByCoordinate \
                    --outFileNamePrefix {self.prefix}"
        cmd2 = f"bamToBed -i star/{self.prefix}Aligned.sortedByCoord.out.bam -bed12 > {self.prefix}.bed"
        cmd3 = f"bedToGenePred {self.prefix}.bed /dev/stdout| genePredToGtf file /dev/stdin  tmp"
        cmd4 = f"sed -i 's/\/dev\/stdin/{self.prefix}_genome/g'  tmp"    # sed 
        cmd5 = '''awk -F '\t' -v OFS='\t' '{if ($3=="transcript") $3="gene"; print $0}' tmp  > tmp.gtf'''  # transcript to gene
        #cmd5 = "awk -F '\t' -v OFS='\t' '{if ($3==\"transcript\") $3=\"gene\"; print $0}' tmp  > tmp.gtf"  # another way
        cmd6 = f"mv tmp.gtf {self.prefix}.gtf"
        cmd7 = f"rm tmp"

        if self.update == False:
            os.chdir("star")
            subprocess.check_call(cmd1, shell = True)
            os.chdir("../")
            subprocess.check_call(cmd2, shell = True)
        
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)
        subprocess.check_call(cmd6, shell = True)
        subprocess.check_call(cmd7, shell = True)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--genome', help='genome', required=True)
    parsers.add_argument('--fasta', help='fasta', required=True)
    parsers.add_argument('--prefix', help='prefix', required=True)
    parsers.add_argument('--star_version', help='star_version', default="old")
    parsers.add_argument('--update', action='store_true')

    args = parsers.parse_args()
    runner = Fasta_to_gtf(args) 
    runner.run()

if __name__ == '__main__':
    main()





