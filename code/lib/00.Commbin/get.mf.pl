use strict;
use warnings;

die"usage:perl $0 <in:group> <out:mf>\n" unless @ARGV ==2;
my($group,$mf)=@ARGV;

open(MF,">$mf");
&write_mf($group,"mf",*MF);
close MF;

sub write_mf{
#=============
    my ($group,$mark,$mf) = @_;
    ($group && -s $group) || die $!;
    open GP,$group || die$!;
    print $mf  "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n";
    while(<GP>){
        my ($sname,$gname) = split/\s+/;
        if($mark eq 'mf'){
        	print $mf "$sname\t\t\t$gname\n";
	}elsif(	$mark eq 'list' ){
		print $mf "$sname\t\t\t$sname\n";
	}else{die "choosing from mf|list for write mf!\n";}
    }
    close GP;
}
