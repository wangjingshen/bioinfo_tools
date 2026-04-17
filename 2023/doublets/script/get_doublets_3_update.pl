#!usr/bin/perl -w
use strict;
use Cwd;
#use File::Basename qw(basename dirname);
my $wd = getcwd;
system"mkdir -p doublets";
my @mtx = `ls ./*/0[56].count/*_downsample.t*`;
my $cmd1;
foreach (@mtx){
chomp;
if(/.\/(\S+)?\/.*\/.*/) {
$cmd1 = << "END";
source  activate celescope1.4.1
Rscript /SGRNJ06/randd/USER/wangjingshen/script/doublets/doublet.r $wd/$1/05.count/$1_filtered_feature_bc_matrix/ $wd/$1/05.count/$1_matrix.tsv   # fix
gzip --force  $wd/$1/05.count/$1_matrix.tsv
cp $wd/$1/05.count/$1_filtered_feature_bc_matrix/genes.tsv $wd/$1/05.count/$1_filtered_feature_bc_matrix/genes.tsv.bak
sed -i "1igene_id\\tgene" $wd/$1/05.count/$1_filtered_feature_bc_matrix/genes.tsv.bak
gunzip -dc $wd/$1/05.count/$1_matrix.tsv.gz| paste $wd/$1/05.count/$1_filtered_feature_bc_matrix/genes.tsv.bak - |awk '{FS="\\t";OFS="\\t";\$2="";\$3="";print \$0}' |sed 's/\\t\\t/\\t/g' > $wd/$1/05.count/$1_matrix.tsv
gzip --force $wd/$1/05.count/$1_matrix.tsv
source activate old
perl /SGRNJ06/randd/USER/wangjingshen/script/doublets/mix_UMI_split_hz.v2.pl \\
  -input  $wd/$1/05.count/$1_matrix.tsv.gz \\
  -name $1 \\
  -out $wd/doublets\n
END
open (RUN,"> $wd/log/$1_doublets.sh")|| die $!;
print RUN "$cmd1";
close RUN;
system "qsub -cwd -V -q randd.q -l vf=1g,p=1 $wd/log/$1_doublets.sh";
}
}
