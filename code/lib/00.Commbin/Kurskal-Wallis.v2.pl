#http://blog.csdn.net/anshiquanshu/article/details/68954536
use strict;
use warnings;

die"usage:perl $0 <in:table> <in:group> <output:prefix> <opts:cutoff> \n" unless @ARGV >=3;
my($table,$group_in,$output)=@ARGV;
my$cutoff=$ARGV[3] if $ARGV[3]; 

my %s2group;
my (@groups,@samples);
open OR,$group_in;
while (<OR>) {
	chomp;
	next if /^#/;
	my@or=split/\s+/;
	$s2group{$or[0]}=$or[1];
	push @groups,$or[1] if !grep{$or[1] eq $_ } @groups;
	push @samples,$or[0];
}
close OR;

my (%s2abun,%out2out,%out2avg,%out2table,%out2box,@tax_p,@pvalue_p);
open OR,"$table";
my $head=<OR>;chomp$head;
my @head=split/\t/,$head;pop@head;
while (my $l=<OR>) {
	chomp$l;my@or=split/\t/,$l;
	next if $or[0] eq 'Other';
	my (%g2abun);
	my $nzeronum=0;
	for (my $i = 1; $i < @or -1; $i++){
		next unless $s2group{$head[$i]};
		if ($or[$i]) {
			push @{$g2abun{$s2group{$head[$i]}}},$or[$i];
			$s2abun{$or[0]}{$head[$i]}=$or[$i];
			$nzeronum++ if $or[$i] ne '0.0';
		}else{
			push @{$g2abun{$s2group{$head[$i]}}},"0";
		}
	}
	next if $cutoff && $nzeronum < $cutoff;
	my (@group,@data,$flag);
	$out2out{$or[0]}="$or[0]";
	$out2avg{$or[0]}="$or[0]";
	foreach my $g (@groups){
		push @group,($g) x @{$g2abun{$g}};
		push @data,@{$g2abun{$g}};
		#print "$g\n@{$g2abun{$g}}\n" if $or[0] eq "Pseudomonadaceae";
		my $median=&median(@{$g2abun{$g}});
		my $avg=&avg(@{$g2abun{$g}});
		$out2out{$or[0]}.="\t$median";
		$median=sprintf("%.2e", $median);
		$avg=sprintf("%.2e", $avg);
		$median=0 if $median eq "0.00e+00";
		$avg=0 if $avg eq "0.00e+00";
		$out2avg{$or[0]} .= "\t$median($avg)";
		if($median == 0 ){$flag .= 0;}else{$flag .= 1;}
	}
	#print "$or[0]\t$flag\n";
	next if $flag !~ /1/;
	my $group="\"".join("\",\"",@group)."\"";
	my $data=join(",",@data);
	my $pvalue=`Rscript --vanilla --slave -e 'group <- c($group);data <- c($data);data85 <- data.frame(group,data);p <- kruskal.test(data85\$data~data85\$group); p\$p.value;'`;
	$pvalue=~s/^\[1\]\s+//;$pvalue=~s/\s+//g;
	#$out2out{$or[0]}.="\t$pvalue";
	if ($pvalue && $pvalue < 0.05) {
		push @tax_p,$or[0];
		push @pvalue_p,$pvalue;
		for (my $i = 1; $i < @or -1; $i++){
			next unless $s2group{$head[$i]};
			$out2box{$or[0]} .= "$or[0]\t$or[$i]\t$s2group{$head[$i]}\t$head[$i]\n";
		}
	}
}
close OR; 

my $pvalues=join(",",@pvalue_p);
#using fdr for correction
my $qvalues=`Rscript --vanilla --slave -e 'p <- c($pvalues);p.adjust(p, method = "fdr", n = length(p));'`;
$qvalues=~s/\[\d+\]//g;
$qvalues=~s/^\s+//;
my @qvalues=split/\s+/,$qvalues;
my (%tax2qvalue,%tax2pvalue);
for(my $i=0;$i < @tax_p;$i++){$tax2qvalue{$tax_p[$i]}=$qvalues[$i]; $tax2pvalue{$tax_p[$i]}=$pvalue_p[$i];}

open OUT,">$output.test.txt";
open OUT2,">$output.test2.txt";
open TABLE,">$output.table.sig.txt";
open BOX,">$output.box.sig.txt";
my @head_groups=map{$_.".(median)";} @groups;
my @head_groups2=map{$_.".(median)(avg)";} @groups;
print OUT "ID\t".join("\t",@head_groups),"\tp.value\tq.value\n";
print OUT2 "ID\t".join("\t",@head_groups2),"\tp.value\tq.value\n";
print TABLE  "ID\t".join("\t",@samples)."\n";
print BOX "ID\tAbun\tGroup\tSample\n";
foreach my $tax (@tax_p){
	next if $tax2qvalue{$tax} && $tax2qvalue{$tax} > 0.05;
	print OUT $out2out{$tax}."\t$tax2pvalue{$tax}\t$tax2qvalue{$tax}\n";
	print OUT2 $out2avg{$tax}."\t$tax2pvalue{$tax}\t$tax2qvalue{$tax}\n";
	print BOX $out2box{$tax};
	print TABLE "$tax";
	foreach my $s (@samples){
		$s2abun{$tax}{$s} ? 
		print TABLE "\t$s2abun{$tax}{$s}":
		print TABLE "\t0";
	}
	print TABLE "\n";
}
close TABLE;
close BOX;
close OUT;
close OUT2;

sub avg{
	my($num,$total);
	foreach my $n (@_){
		$num++;
		$total += $n;
	}
	return $total/$num;
}

sub median{
    my @list = sort{$a<=>$b} @_;
    my $count = @list;
    if( $count == 0 )
    {
        return undef;
    }   
    if(($count%2)==1){
        return $list[int(($count-1)/2)];
    }
    elsif(($count%2)==0){
        return ($list[int(($count-1)/2)]+$list[int(($count)/2)])/2;
    }
}
