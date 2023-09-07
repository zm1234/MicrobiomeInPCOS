use strict;
use warnings;
use PerlIO::gzip;
use Cwd qw(abs_path);

die"usage:perl $0 <in:data.list> <outdir>\n" unless @ARGV == 2;
my($in,$outdir)=@ARGV;
-s $outdir || `mkdir -p $outdir`;
$outdir=abs_path($outdir);

my (%fqs,@samples);
open IN,$in;
while(<IN>){
	chomp;
	my($s,$fq)=(split/\s+/)[0,1];
	my $flag=$1 if $fq=~/.*_[12].fq.gz$/;
	push @{$fqs{$s}{$flag}},$fq;
	push @samples,$s;
}
close IN;

open OUT,">$outdir/raw.data.list";
foreach my $s (@samples){
	print OUT "$s\t$outdir/$s\_raw_1.fq.gz,$outdir/$s\_raw_2.fq.gz\n";
	open O,">:gzip","$outdir/$s\_raw_1.fq.gz";
	foreach my $file (@{$fqs{$s}{1}}){
		$file=~/\.gz$/ ? open FI,"gzip -dc $file|"  : open FI,$file;
		while(my $l=<FI>){print O $l;}
		close FI;
	}
	close O;

	open O,">:gzip","$outdir/$s\_raw_2.fq.gz";
	foreach my $file (@{$fqs{$s}{1}}){
		$file=~/\.gz$/ ? open FI,"gzip -dc $file|"  : open FI,$file;
		while(my $l=<FI>){print O $l;}
		close FI;
	}
	close O;
}
close OUT;
