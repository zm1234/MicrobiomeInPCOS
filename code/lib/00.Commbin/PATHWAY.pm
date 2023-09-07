package  PATHWAY;
use strict qw(subs refs);
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw( get_pathway );

sub get_pathway{
	my($file,$key_words)=@_;
	open FILE,$file;
	my %key2path;
	while (my $l=<FILE>) {
		chomp$l;
		if ($l=~/(\S+)\s+==\s+(.*)/) {
			my $key=$1;
			my $path=$2;
			$key2path{$key}=$path;
		}elsif($l=~/(\S+)\s+=\s+(.*)/){
			my $key=$1;
			my $path=$2;
			foreach my $k (keys %key2path){
				$path=~s/$k/$key2path{$k}\//g;
			}
			$key2path{$key}=$path;
		}
	}
	close FILE;

	my @result;
	foreach my $key (@$key_words){
		$key2path{$key} ? 
		push @result,$key2path{$key} :
		push @result,"";
	}
	return @result;
}

1;
