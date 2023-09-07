use strict;
use warnings;

die"usage:perl $0 <in:otu table> <in:mf> <out:otu table> <out:tags.stat>\n" unless @ARGV == 4;
my($in,$mf,$out_table,$out_stat)=@ARGV;

open(OR,$in);
open(OUT,">$out_table");
my $line1=<OR>;
print OUT $line1;
chomp(my$head=<OR>);
my@head=split/\t/,$head;

print OUT "$head\n";
my %sample2info;
while (<OR>) {
	chomp;
	my@or=split/\t/;
	my $detail=pop@or;
	my $otu_id=shift@or;
	my $total_tags;
	map{ $_=~s/\.0$//;$total_tags+=$_;} @or;

	for (my $i = 0; $i < @or; $i++) {
		$sample2info{$head[$i+1]}{Total_Tags} += $or[$i];
	}

	if($total_tags > 1){
		print OUT "$otu_id\t".join("\t",@or)."\t$detail\n";
		for (my $i = 0; $i < @or; $i++) {
			$sample2info{$head[$i+1]}{OTU_num} ++ if $or[$i] ;
		}
		if($detail eq 'Unassigned') {
			for (my $i = 0; $i < @or ; $i++){
				$sample2info{$head[$i+1]}{Unassigned} += $or[$i];
			}
		}else{
			for (my $i = 0; $i < @or ; $i++){
				$sample2info{$head[$i+1]}{Taxon_Tags} += $or[$i];
			}
		}
	}else{
		for (my $i = 0; $i < @or; $i++) {
			$sample2info{$head[$i+1]}{uniq_tags} += $or[$i];
		}
	} 	
}
close OR;
close OUT;

open(OUT,">$out_stat");
print OUT "Sample_Name\tTotal_tag\tTaxon_Tag\tUnassigned_Tag\tUnique_Tag\tOTU_num\n";
open(OR,$mf);
while (<OR>) {
	chomp;
	next if /^#/;
	my@or=split/\s+/;
	next if ! $sample2info{$or[0]};
	$sample2info{$or[0]}{uniq_tags} ||= 0;
	$sample2info{$or[0]}{Unassigned} ||= 0;
	my@comma_add=($sample2info{$or[0]}{Total_Tags},$sample2info{$or[0]}{Taxon_Tags},$sample2info{$or[0]}{Unassigned},$sample2info{$or[0]}{uniq_tags},$sample2info{$or[0]}{OTU_num});
	&comma_add(@comma_add);
	print OUT "$or[0]\t".join("\t",@comma_add)."\n";
}
close OR;
close OUT;

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

