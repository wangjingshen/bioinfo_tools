import os
import argparse
import subprocess

class Cp_lims_mtx:
    def __init__(self, args):
        self.expM = args.expM
        self.name = args.name
        self.mode = args.mode
        #self.new_path = args.new_path
        #if not os.path.exists(f"{self.new_path}"):
        #    os.system(f"mkdir -p {self.new_path}")

    def run(self):
        if self.mode == "cloud":
            cmd1 = f"mkdir -p {self.name}"
            cmd2 = f"cp {self.expM} {self.name}"
            cmd3 = f"tar -xvf {self.name}/{self.name}_matrix_10X.tar -C {self.name}"

        if self.mode == "local":
            cmd1 = f"cp -r {self.expM}/ {self.name}"
            cmd2 = f"tar -cf {self.name}_matrix_10X.tar --directory {self.name} ."
            cmd3 = f"mv {self.name}_matrix_10X.tar {self.name}"

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)



def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--expM', help='expM', required=True)
    parsers.add_argument('--name', help='name', required=True)
    parsers.add_argument('--mode', help='mode', required=True)

    args = parsers.parse_args()
    runner = Cp_lims_mtx(args) 
    runner.run()

if __name__ == '__main__':
    main()