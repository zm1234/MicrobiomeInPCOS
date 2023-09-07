use strict;
use warnings;
use PerlIO::gzip;

die "usage:perl $0 <in:fq.gz> <cutoff> <out> \n" unless @ARGV == 3;
my($in,$cut,$out)=@ARGV;

$in=~/.gz$/ ? open IN,"gunzip -dc $in |" : open IN,$in || die $!;
#open IN,"gzip -dc $in|" : $in=~/\.gz$/ ? open IN, "<:gzip", "$in" : open IN,$in || die $!;
$out=~/\.gz$/ ? open OUT,">:gzip", "$out" : open OUT,">$out" || die $!;
while (my $head=<IN>){
	my $seq=<IN>;chomp$seq;
	<IN>;
	my $qua=<IN>;chomp$qua;
	my $len=length($seq)-$cut;
	my $seq_new=substr($seq,0,$len);
	my $qua_new=substr($qua,0,$len);
	print OUT "$head$seq_new\n+\n$qua_new\n";
}
close IN;
close OUT;
