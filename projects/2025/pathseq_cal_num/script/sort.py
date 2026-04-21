import subprocess
import sys

def sort(spname):
    cmd1 = f'head -n 1 {spname}.tsv > tmp1'
    cmd2 = f'tail -n +2 {spname}.tsv | sort -k 2 -n -r > tmp2'
    cmd3 = f'cat tmp1 tmp2 > {spname}_sort.tsv'
    cmd4 = f'rm tmp1 tmp2'

    subprocess.check_call(cmd1, shell = True)
    subprocess.check_call(cmd2, shell = True)
    subprocess.check_call(cmd3, shell = True)
    subprocess.check_call(cmd4, shell = True)

def main():
    spname = sys.argv[1].split(',')
    results = [sort(value) for value in spname]


if __name__ == '__main__':
    main()