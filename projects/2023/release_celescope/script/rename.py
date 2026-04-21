import os
import argparse
import subprocess

class Rename:
    def __init__(self, args):
        self.old_name = args.old_name
        self.new_name = args.new_name

    def run(self):
        cmd1 = f"rename -v {self.old_name} {self.new_name} {self.old_name} 1>rename.o 2>rename.e"
        cmd2 = f"rename -v {self.old_name} {self.new_name} {self.new_name}/*/*{self.old_name}* 1>>rename.o 2>>rename.e"
        cmd3 = f"rename -v {self.old_name} {self.new_name} {self.new_name}/*{self.old_name}* 1>>rename.o 2>>rename.e"
        cmd4 = f"sed -i 's/{self.old_name}/{self.new_name}/g' `grep -rl {self.old_name} {self.new_name}/*/*`"
        cmd5 = f"sed -i 's/{self.old_name}/{self.new_name}/g' {self.new_name}/.data.json"
        cmd6 = f"sed -i 's/{self.old_name}/{self.new_name}/g' {self.new_name}/.metrics.json"
        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)
        subprocess.check_call(cmd6, shell = True)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--old_name', help='genome', required=True)
    parsers.add_argument('--new_name', help='new_name', required=True)

    args = parsers.parse_args()
    runner = Rename(args) 
    runner.run()

if __name__ == '__main__':
    main()





