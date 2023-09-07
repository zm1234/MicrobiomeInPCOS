use strict;
use warnings;

die"usage:perl $0 <indir> <inmf> <outfile>\n" unless @ARGV ==3;

my($indir,$mf,$out)=@ARGV;

my %sample2info;
my $chimedir;
my $flashdir;
foreach my $file (split/\s+/,`ls $indir/*/*.fqcheck`){
	next if ! $file;
	$file=~s/\/\//\//g;
	my$sample=$1 if $file=~/.*\/(.*).fqcheck/;
	$sample2info{$sample}{lastdata}=(split/\s+/,`grep 'Readnum' $file`)[1];
	$sample2info{$sample}{datasize}=(split/\s+/,`grep 'DataSize' $file`)[1];
	$sample2info{$sample}{q20}=(split/\s+/,`grep 'Q20' $file`)[1];
	$sample2info{$sample}{q30}=(split/\s+/,`grep 'Q30' $file`)[1];
	$sample2info{$sample}{gc}=(split/\s+/,`grep 'GC' $file`)[1];
	$sample2info{$sample}{avglen}=$sample2info{$sample}{lastdata} ? sprintf("%d",$sample2info{$sample}{datasize}/$sample2info{$sample}{lastdata}) : "0";
	$chimedir="$indir/$sample/03.uchime";
	$flashdir="$indir/$sample/01.combine";
}

$chimedir=~s/\/\//\//g;
$flashdir=~s/\/\//\//g;

if(-s $flashdir){
	foreach my $file (split/\s+/,`ls $indir/*/01.combine/*.flash.log`){
		next if ! $file;
		$file=~s/\/\//\//g;
		my$sample=$1 if $file=~/.*\/(.*).flash.log/;
		my $total=(split/\s+/,`grep 'Total reads' $file`)[3];
		my $combine=(split/\s+/,`grep 'Combined reads' $file`)[3];
		$sample2info{$sample}{'total'}=$total;
		$sample2info{$sample}{'combine'}=$combine;
	}

	foreach my $file (split/\s+/,`ls $indir/*/02.qc/split_library_log.txt`){
		next if ! $file;
		$file=~s/\/\//\//g;
		my$sample=$1 if $file=~/.*\/(.*)\/02.qc\/split_library_log.txt/;
		my $qcdata=(split/\s+/,`grep 'Total number seqs written' $file`)[4];
		$sample2info{$sample}{'qcdata'}=$qcdata;
	}
}else{
	foreach my $file (split/\s+/,`ls $indir/*/02.qc/split_library_log.txt`){
		next if ! $file;
		$file=~s/\/\//\//g;
		my$sample=$1 if $file=~/.*\/(.*)\/02.qc\/split_library_log.txt/;
		my $total=(split/\s+/,`grep 'Total number of input sequences' $file`)[5];
		my $qcdata=(split/\s+/,`grep 'Total number seqs written' $file`)[4];
		$sample2info{$sample}{'total'}=$total;
		$sample2info{$sample}{'qcdata'}=$qcdata;
	}
}

if (-s $chimedir) {
	open(OR,$mf);
	open(OUT,">$out");
	-s $flashdir ? 
	print OUT "Sample Name\tRaw PE(#)\tRaw Tags(#)\tClean Tags(#)\tEffective Tags(#)\tBase(nt)\tAvglen(nt)\tQ20\tQ30\tGC%\tEffective%\n":
	print OUT "Sample Name\tRaw Reads(#)\tClean Reads(#)\tEffective Reads(#)\tBase(nt)\tAvglen(nt)\tQ20\tQ30\tGC%\tEffective%\n";
	while (<OR>) {
		next if /^#/;
		my($sample)=(split/\s+/)[0];
		$sample2info{$sample}{'lastdata'} ||= 0;
		$sample2info{$sample}{'total'} ||= 0;
		$sample2info{$sample}{'combine'} ||= 0;
		$sample2info{$sample}{'qcdata'} ||= 0;
		$sample2info{$sample}{'datasize'} ||= 0;
		
		my $effective=$sample2info{$sample}{'lastdata'} ? 
		sprintf("%.2f",$sample2info{$sample}{'lastdata'}/$sample2info{$sample}{'total'}*100).'%' : 0;
		my @comma_add= 
		-s $flashdir ? 
		($sample2info{$sample}{'total'},$sample2info{$sample}{'combine'},$sample2info{$sample}{'qcdata'},$sample2info{$sample}{'lastdata'},$sample2info{$sample}{'datasize'}) :
		($sample2info{$sample}{'total'},$sample2info{$sample}{'qcdata'},$sample2info{$sample}{'lastdata'},$sample2info{$sample}{'datasize'});
		&comma_add(@comma_add);
		print OUT "$sample\t".join("\t",@comma_add)."\t$sample2info{$sample}{'avglen'}\t$sample2info{$sample}{'q20'}\t$sample2info{$sample}{'q30'}\t$sample2info{$sample}{'gc'}\t$effective\n";
	}
	close OR;
	close OUT;
}else{
	open(OR,$mf);
	open(OUT,">$out");
	-s $flashdir ?
	print OUT "Sample Name\tRaw PE(#)\tRaw Tags(#)\tClean Tags(#)\tBase(nt)\tAvglen(nt)\tQ20\tQ30\tGC%\tEffective%\n":
	print OUT "Sample Name\tRaw Reads(#)\tClean Reads(#)\tBase(nt)\tAvglen(nt)\tQ20\tQ30\tGC%\tEffective%\n";
	while (<OR>) {
		next if /^#/;
		my($sample)=(split/\s+/)[0];
		my $effective=sprintf("%.2f",$sample2info{$sample}{'lastdata'}/$sample2info{$sample}{'total'}*100).'%';
		my @comma_add=
		-s $flashdir ?
		($sample2info{$sample}{'total'},$sample2info{$sample}{'combine'},$sample2info{$sample}{'lastdata'},$sample2info{$sample}{'datasize'}):
		($sample2info{$sample}{'total'},$sample2info{$sample}{'lastdata'},$sample2info{$sample}{'datasize'});
		&comma_add(@comma_add);
		print OUT "$sample\t".join("\t",@comma_add)."\t$sample2info{$sample}{'avglen'}\t$sample2info{$sample}{'q20'}\t$sample2info{$sample}{'q30'}\t$sample2info{$sample}{'gc'}\t$effective\n";
	}
	close OR;
	close OUT;
}
##########################sub routines####################################
sub comma_add{
    foreach(@_){
        $_ || next;
        $_ = /(\d+)\.(\d+)/ ? comma($1) . '.' . comma($2,1) : comma($_);
    }
}
sub comma{
    my ($c,$rev) = @_;
    (length($c) > 3) || return($c);
    $rev || ($c = reverse $c);
    $c =~ s/(...)/$1,/g;
    $rev || ($c = reverse $c);
    $c =~ s/^,//;
    $c;
}
