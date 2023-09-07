#Junru Chen, chenjr@geneworks.cn
##2016-10-21
use strict;
use warnings;

die "usage:perl $0 <in:> <out>\n" unless @ARGV == 2;
my($in,$out)=@ARGV;

open(OR,"$in");
my$head=<OR>;chomp$head;
my@head=split/\t/,$head;
shift@head;shift@head;shift@head;
my%hash;
while(<OR>){
	chomp;
	my@or=split/\t/;
	shift@or;
	my$sequence_num=shift@or;
	my$iter=shift@or;
	for(my$i=0;$i < @or;$i++){
		$hash{$sequence_num}{$head[$i]}{'total'} += $or[$i];
		$hash{$sequence_num}{$head[$i]}{'num'} ++;	
	}
}
close OR;

open(OUT,">$out");
print OUT "sequence_num\tsample\tindex\n";	
foreach my $num (sort {$a <=> $b} keys %hash){
	foreach my $s (@head){
		print OUT "$num\t$s\t".$hash{$num}{$s}{'total'}/$hash{$num}{$s}{'num'}."\n";
	}
}

