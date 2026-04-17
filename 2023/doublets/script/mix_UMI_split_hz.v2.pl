#!/usr/bin/perl -w
#use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);
use List::Util qw(max min sum maxstr minstr shuffle);
my $programe_dir=basename($0);
my $path=dirname($0);
my $ver    = "1.0";
my $Writer = "hulongfei\@singleronbio.com";
my $Data   = "2018/06/30";
my $BEGIN=time();
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($infile,$outdir,$name,$umi_cutoff);
GetOptions(
			"h|?" =>\&help,
			"input:s"=>\$infile,
			"name:s"=>\$name,
			"umicf:s"=>\$umi_cutoff,
			"out:s"=>\$outdir,
			) || &help;
&help unless ($infile && $name);
my $dir=getcwd;
$outdir||=$dir;
$umi_cutoff||=0;    # js mod $umi_cutoff||=200
`mkdir $outdir` unless (-e $outdir);
sub help
{
	print <<"	Usage End.";
	Description:
	Writer  : $Writer
	Data    : $Data
	Version : $ver
	function: ......
	Usage:
		-input		input file			must be given
	
		-name		name

		-umicf		UMI cutoff			[200]

		-out		outdir				must be given

		-h		Help				document
	Usage End.
	exit;
}
my %split;
my $bn;
open (G,"> $outdir/$name.mix_Species.txt") || die $!;
if($infile =~ /\.gz$/){
	open (IN,"gzip -dc $infile|") || die $!;
}else{
	open (IN,"$infile")|| die $!;
}
my $head=<IN>;
chomp($head);
$head=~s/\s+$//g;
if($head=~/^geneID/){
	$head=~s/^geneID\s+//g;
}
my @li=split(/\s+/,$head);
my $flag="mix";
while(<IN>){
	chomp;
	s/\s+$//g;
	my @li1=split(/\s+/,$_);
	if($li1[0]=~/^ENSG/){
		$flag="human";
	}elsif($li1[0]=~/^ENSMUSG/){
		$flag="mouse";
	}else{
		warn "$li1[0] is error\n";
		$flag="error";
	}
	for($bn=1;$bn<@li1;$bn++){
		my $i=$bn-1;
		$split{$li[$i]}{$flag}+=$li1[$bn];
	}
}
close IN;
print G "cell_id\thuman\tmouse\tSpecies\n";
my %number;
my $tol=0;
my %com_hm;
my $un_cell=0;
foreach my $bc (@li){
	my $homo;
	my $mus;
	if(! exists $split{$bc}){
		print "$bc is miss\n";
	}
	if(exists $split{$bc}{'human'}){
		$homo=$split{$bc}{'human'};
	}else{
		$homo=0;
	}
	if(exists $split{$bc}{'mouse'}){
		$mus=$split{$bc}{'mouse'};
	}else{
		$mus=0;
	}
	my $type="indeterminate";
	my $total = $homo+$mus;
=pod
	if($total>=2000){
		if($homo>=2000 and $mus < 1000){
			$type="human";
		}
		if($homo<1000 and $mus>= 2000){
			$type="mouse";
		}
	}else{
		$type="indeterminate";
	}
	if($total>=5000 and $mus >= 1000  and $homo >=1000){
		$type="multiplet";
	}
=cut
	#if($total>=50000){    # js mod, save >=50000
	#	next;
	#}
	if($total <$umi_cutoff){
		$un_cell++;
		next;
	}
	my $rato=$homo/$total;
	if($rato > 0.8){
		$type="human";
	}elsif($rato > 0.2 and $rato <= 0.8){
		$type="multiplet";
	}else{
		$type="mouse";
	}
	$number{$type}++;
	$com_hm{$type}+=$total;
	if(exists $com_hm{$type}{'mouse'}){
		$com_hm{$type}{'mouse'}+=$mus;
	}else{
		$com_hm{$type}{'mouse'}=$mus;
	}
	if(exists $com_hm{$type}{'human'}){
		$com_hm{$type}{'human'}+=$homo;
	}else{
		$com_hm{$type}{'human'}=$homo;
	}
	print G "$bc\t$homo\t$mus\t$type\n";
	$tol++;
}
close G;
my $dr=$number{multiplet}*100/($tol);
open (REPORT,"> $outdir/$name.report.xls")|| die $!;
foreach my $ty (sort keys %number){
	print  REPORT  "number:$ty\t$number{$ty}\n";
}
print REPORT  "indeterminate:$un_cell\n";
print REPORT  "doublets:$name\t$dr\n";
my $m2h=0;
my $h2m=0;
if((exists $com_hm{"human"}{"mouse"})  and (exists $com_hm{"human"})){
	$m2h=$com_hm{"human"}{"mouse"}*100/$com_hm{"human"};
}
if((exists $com_hm{"mouse"}{"human"})  and (exists $com_hm{"mouse"})){
	$h2m=$com_hm{"mouse"}{"human"}*100/$com_hm{"mouse"};
}
print REPORT "mouse to human: $m2h\n";
print REPORT "human to mouse: $h2m\n";
#$outdir/mix_Species.txt
open(R,"> $outdir/mix_plot2.v2.r") || die $!;
print R <<"	end.";
#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(argparser)
library(MASS)
library(Cairo)
library(scales)
stat <- read.table("$outdir/$name.mix_Species.txt",sep='\\t',header=T)
stat\$human <- stat\$human/1000
stat\$mouse <- stat\$mouse/1000
levels(stat\$Species)[levels(stat\$Species)=="human"] <- "human($number{human})"
levels(stat\$Species)[levels(stat\$Species)=="mouse"] <- "mouse($number{mouse})"
levels(stat\$Species)[levels(stat\$Species)=="multiplet"] <- "multiplet($number{multiplet})"
p<-ggplot(data=stat,aes(x=human,y=mouse,colour=Species))+ geom_point(size=1)+scale_colour_manual(values=c("red","blue","purple"))+theme_bw()+theme(legend.position=c(1,1),legend.justification=c(1,1),legend.background=element_blank(),legend.key=element_blank(),legend.title=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.title.y=element_text(size =12,face ="bold"),axis.text.x = element_text(size =12,face ="bold"),axis.text.y=element_text(size=12,face ="bold"),axis.title.x=element_text(size=12,face="bold"),plot.title=element_text(size=12,face="bold",hjust = 0.5))+ylab("Number of mouse transcripts\\n(10^3)")+xlab("(10^3)\\nNumber of human transcripts")+ggtitle(paste("$name","",sep=''))
pdf(paste("$outdir",'/',"$name",'_homo_mus.pdf',sep=''))
p
dev.off()
CairoSVG(paste("$outdir",'/',"$name",'_homo_mus.svg',sep=''))
p
dev.off()
	end.
close R;
`Rscript $outdir/mix_plot2.v2.r`;

