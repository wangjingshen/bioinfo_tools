import os
import argparse
import subprocess

class Release_celescope:
    def __init__(self, args):
        self.celescope_path = args.celescope_path
        self.new_path = args.new_path
        self.old_name = args.old_name
        self.new_name = args.new_name

        if not os.path.exists(self.new_path):
            os.system(f'mkdir -p {self.new_path}')

    def ln(self):
        cmd1 = f"ln -s {self.celescope_path} {self.new_path}/{self.new_name}"
        subprocess.check_call(cmd1, shell = True)

    def cp(self):
        cmd1 = f"cp -r {self.celescope_path} {self.new_path}/{self.new_name}"
        subprocess.check_call(cmd1, shell = True)

    def rename(self):
        cmd1 = f"rename -v {self.old_name} {self.new_name} {self.new_path}/{self.new_name}/*/*"   # 1>>rename.o 2>>rename.e"
        cmd2 = f"rename -v {self.old_name} {self.new_name} {self.new_path}/{self.new_name}/*"     # 1>>rename.o 2>>rename.e"
        cmd3 = f"sed -i 's/{self.old_name}/{self.new_name}/g' `grep -rl {self.old_name} {self.new_path}/{self.new_name}/*`"
        cmd4 = f"sed -i 's/{self.old_name}/{self.new_name}/g' {self.new_path}/{self.new_name}/.data.json"
        cmd5 = f"sed -i 's/{self.old_name}/{self.new_name}/g' {self.new_path}/{self.new_name}/.metrics.json"

        if( self.new_name != self.old_name):
            subprocess.check_call(cmd1, shell = True)
            subprocess.check_call(cmd2, shell = True)
            subprocess.check_call(cmd3, shell = True)
            subprocess.check_call(cmd4, shell = True)
            subprocess.check_call(cmd5, shell = True)
    
    def gzip(self):
        cmd1 = f"gzip {self.new_path}/{self.new_name}/01.barcode/{self.new_name}_2.fq"
        cmd2 = f"gzip {self.new_path}/{self.new_name}/02.cutadapt/{self.new_name}_clean_2.fq"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)

    def rm(self):
        cmd1 = f"rm {self.new_path}/{self.new_name}/01.barcode/{self.new_name}_2.fq"
        cmd2 = f"rm {self.new_path}/{self.new_name}/02.cutadapt/{self.new_name}_clean_2.fq"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)

    def run(self):
        if self.new_name == self.old_name:
            self.ln()
        else:
            self.cp()
            self.rename()
            #self.gzip()
            self.rm()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--celescope_path', help='celescope path', required=True)
    parsers.add_argument('--new_path', help='new path', required=True)
    parsers.add_argument('--old_name', help='old name', required=True)
    parsers.add_argument('--new_name', help='new_name', required=True)

    args = parsers.parse_args()
    runner = Release_celescope(args) 
    runner.run()

if __name__ == '__main__':
    main()





