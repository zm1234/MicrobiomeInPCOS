#!/usr/bin/perl -w

#====================================================================================================================
###################################################################################################


use strict;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use lib "$Bin/../lib/00.Commbin";
use PATHWAY;
use Getopt::Long;

#set default paramers
my %opt = (
    "type","pe","db","greengene","otu","uparse","region","16S_V4","step2",12345,
    #options for step1
    #"flash", " -m 10 -f 300 -x 0.1 -p 33 -r 148 -M 71 ",
    "qc","-q 19 --barcode_type not-barcoded --store_demultiplexed_fastq",
    #"primer",20,
    "primer",13, #change by songde
    "adapt","AGATCGGAAGAGC", #add songdq in 20190906
    "ng_QC","-N 0.1 -q 33 -L 5 -p 0.5",
    #options for step2
    "identity",0.97,"top",10,"heatmap",30,
    #other options
    "outdir","./","step",12345,"lda",2
);

#get options
GetOptions(
    \%opt,
    "data_list:s","type:s","mf:s","db:s","otu:s","region:s","cut:s","primer:s","group:s","venn:s","venn_group:s","step2:s",
    #options for step1
    "flash:s","qc:s","ng_QC:s","adapt:s",
    #options for step2
    "identity:s","top:s","heatmap:s",
    #other options
    "outdir:s","shdir:s","step:s","notrun","help:s",
);

#get software pathway
(-s "$Bin/Pathway_cfg.txt") || die"error: can't find config at $Bin, $!\n";
my($kur_wal,$kw_alpha,$otu2lefse,$lefse,$read_cut,$primer_cut,$Rscript,$fqcheck,$flash,$split_libraries_fastq,$usearch,$step1_lib,$ng_QC,$gold,$step2_lib,$avg_group,$sort_max,$fa_16s,$tax_16s,$qiime,$core_set,$usearch11,$usearch5,$uparse64,$step3_lib,$step4_lib,$silva)=get_pathway("$Bin/Pathway_cfg.txt",[qw(KUR_WAL KW_ALPHA OTU2LEFSE LEFSE READ_CUT PRIMER_CUT Rscript FQCHECK FLASH QIIME_QC USEARCH STEP1_LIB NG_QC GOLD STEP2_LIB AVG_GROUP SORT_MAX 16S_FA 16S_TAX QIIME CORE_SET USEARCH11 USEARCH5 UPARSE64 STEP3_LIB STEP4_LIB SILVA )]);
my $qsub="perl $Bin/../lib/00.Commbin/super.qsub.pl ";
my $fqstat="perl $step1_lib/Fastq.stats.pl";
my $fqsort="perl $step1_lib/sort.fastq.pl";
my $get_result="perl $Bin/../lib/00.Commbin/get.result.pl ";
my $get_report="perl $Bin/../lib/00.Commbin/get.report.pl ";

#====================================================================================================================
###################################################################################################
#options set
die `pod2text $0\n` unless $opt{data_list} && -s $opt{data_list};
$opt{data_list}=abs_path($opt{data_list});
$opt{shdir} = "$opt{outdir}/Shell" if !$opt{shdir};
die "Options set wrong:type must set from pe|se!\n" unless $opt{type} && ($opt{type} eq 'pe' || $opt{type} eq 'se') ;
$flash .= " $opt{flash} " if $opt{flash};
$split_libraries_fastq .= " $opt{qc} " if $opt{qc};
$ng_QC .= " $opt{ng_QC} " if $opt{ng_QC};
die "Database is error!\n\n" unless $opt{db}=~/^greengene/i or $opt{db}=~/^silva/i;
die "OTU method is error!\n\n" unless $opt{otu}=~/^uparse/i or $opt{otu}=~/^usearch/i;
my $assign_taxonomy="$qiime/assign_taxonomy.py ";
$assign_taxonomy .= " --reference_seqs_fp $fa_16s --id_to_taxonomy_fp $tax_16s " if $opt{db}=~/^greengene/i;
$assign_taxonomy .= " --reference_seqs_fp $silva/SILVA_128_SSURef_Nr99_tax_silva_trunc.16S.fasta  --id_to_taxonomy_fp $silva/SILVA_128_SSURef_Nr99_tax_silva_trunc.16S.tax.txt " if $opt{db}=~/^silva/i;
$core_set="$silva/core_alignment_SILVA128.fna" if $opt{db}=~/^silva/i;

$opt{type}=lc($opt{type});
my($r1_cut,$r2_cut)=(split/,/,$opt{cut})[0,1] if $opt{cut};


#====================================================================================================================
###################################################################################################
# judge sample to group 
my $sample_num = judge_sample_group($opt{data_list},$opt{mf});
my ($top_heatmap,$rare,$rank,$pcoa);
if ($sample_num >= 35){
    $top_heatmap="$step2_lib/top_heatmap.large.r";
    $rare="$step3_lib/alpha.diversity.large.R";
    $rank="$step3_lib/Rank.Abundance.large.R";
    $pcoa="$step4_lib/DendoPCoA.large.R";
}else{
    $top_heatmap="$step2_lib/top_heatmap.r";
    $rare="$step3_lib/alpha.diversity.large.R";
    $rank="$step3_lib/Rank.Abundance.R";
    $pcoa="$step4_lib/DendoPCoA.R";
}
#make directorys
my $detail = "$opt{shdir}/detail";
for($opt{outdir},$opt{shdir},$detail){
    (-s $_) || mkdir($_);
    $_ = abs_path($_);
}
my @dir;
for ("01.Combine_QC","02.OTU_analysis","03.Alpha_diversity","04.Beta_diversity","05.Stat_test"){
    push @dir,"$opt{outdir}/$_";
}
#make mf
open(MF,">$opt{shdir}/samples.mf");
$opt{mf} ? 
&write_mf($opt{mf},"mf",*MF) :
&write_mf($opt{data_list},"list",*MF);
$opt{mf}="$opt{shdir}/samples.mf";
close MF;
#====================================================================================================================
###################################################################################################
#main scripts
my $main_shell = "16S_pipeline.pbs";
open SH,">$opt{shdir}/$main_shell" || die$!;
print SH
"#!/bin/bash
#\$ -cwd
#\$ -S /bin/bash
#\$ -N main.16s
#\$ -q NGS\n";
my $step=1;


##==  1) combine and qc for samples  ==##
my $cleanlist;
if($opt{step}=~/1/){
    (-d $dir[0]) || mkdir $dir[0];
    #outdir, shell dir, input data_list
    my ($shell1,$shell2);
    ($cleanlist,$shell1,$shell2)=&write_step1("$dir[0]","$detail/01.Combine_QC","$opt{data_list}","$opt{mf}");
    write_main(*SH,"$step) run Data Control","$detail/01.Combine_QC/",["$shell1.sh","$shell2.sh"],[" --prefix combine_qc --dir $shell1 "," --prefix qc_stat --dir $shell2 "]);
    $step++;
}
#$cleanlist ||= "$dir[0]/cleandata.list";
$cleanlist= "$dir[0]/cleandata.list" if (-s "$dir[0]/cleandata.list"); 
$cleanlist ||= $opt{data_list}; 

##== 2) OTU analysis  ==##
my($otu_table,$otu_biom,$tree,$summarize);
if($opt{step}=~/2/ ){
    (-d $dir[1]) || mkdir $dir[1];
    my ($shell1,$shell2);
    ($otu_table,$otu_biom,$tree,$summarize,$shell1,$shell2)=&write_step2($dir[1],"$detail/02.OTU_analysis",$cleanlist,$opt{mf});
    write_main(*SH,"$step) run OTU analysis","$detail/02.OTU_analysis/",["$shell1.sh","$shell2.sh"],[" --prefix make_otu --dir $shell1 --parallel 1 --ppn 12 "," --prefix otu_stat --dir $shell2 "]);
    $step++;
}
$otu_table ||= "$dir[1]/01.Make_OTU/04.otu_table_with_taxonomy.txt";
$otu_biom ||=  "$dir[1]/01.Make_OTU/04.otu_table.even.biom";
$tree ||=  "$dir[1]/01.Make_OTU/rep_set.tre";
$summarize ||= "$dir[1]/01.Make_OTU/03.summarize.sorted.otu.table.txt";


##== 3) Alpha diveristy  ==##
if($opt{step}=~/3/ ){
    (-d $dir[2]) || mkdir $dir[2];
    my $shell=&write_step3($otu_table,$otu_biom,$tree,$summarize,$opt{mf},"$detail/03.Alpha_diversity",$dir[2]);
    write_main(*SH,"$step) run alpha diveristy","$detail/03.Alpha_diversity/",["$shell.sh"],[" --prefix alpha --dir $shell --parallel 1 "]);
    $step++;
}

##== 4) Diversity Analysis  ==##
if($opt{step}=~/4/ ){
    (-d $dir[3]) || mkdir $dir[3];
    my $shell=write_step4($otu_biom,$tree,$opt{mf},$dir[3],"$detail/04.Beta_diveristy");
    write_main(*SH,"$step) run beta diveristy","$detail/04.Beta_diveristy/",["$shell.sh"],[" --prefix beta --dir $shell --parallel 1 "]);
    $step++;
}


##== 5) Get report result ==##
if($opt{step}=~/5/){
    write_main(*SH,"$step) get result&report","$detail",["$get_result $opt{outdir} $opt{mf} $detail/05.result_report/ $opt{outdir}/\n$get_report $opt{outdir}/result $opt{outdir}"]);
}

close SH;

$opt{notrun} || system"cd $opt{shdir}; qsub $main_shell";


#====================================================================================================================
###################################################################################################
###################################################################################################
##=====================================##
##            SUB FUNCTIONs            ##
##=====================================##
#=============
sub judge_sample_group{
   my($sample_list,$mf)=@_;
   my $num;
   if ($mf and -s $mf){
	my ($num1,$num2,%samples_list,%mf_samples);
	open DD, $sample_list ||die $!;
	while (<DD>){
	    next if (/^#/);
	    my @pp=split/\s+/;
	    $samples_list{$pp[0]}=1;
	    if($pp[1]=~/,/){
		my @files_fqs=split/,/,$pp[1];
	    	die "$pp[0]:$pp[1] is not exists!\n\n" unless -s $files_fqs[0] and  -s $files_fqs[1];
	   }else{
		die "$pp[0]:$pp[1] is not exists!\n\n" unless -s $pp[1];} #add by songdq, cleandata.list 
	    $num1++;
	}
	close DD;
	open AA, $mf||die $!;
	while (<AA>){
	   next if (/^#/);
	   my @pp=split/\s+/;
	   $mf_samples{$pp[0]}=1;
	   $num2++;
	}
	close AA;
	my $flag=1;
	foreach my $i(keys %samples_list){
		if(! exists $mf_samples{$i}){
			$flag=0;
			last;
		}
	}
	foreach my $i(keys %mf_samples){
		if(! exists $samples_list{$i}){
			$flag=0;
			last;
		}
	}
	if ($flag && $num2 == $num1){
	    $num=$num1
	}else{
	   die "Please check samples name of data list file and group file!\n\n";
	}
   }else{
	open DD, $sample_list ||die $!;
	while (<DD>){
	    next if (/^#/);
	    my @pp=split/\s+/;
	    $num++;
	}
	close DD;
   }
	return $num;
}
#=============

#=============
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
    }elsif( $mark eq 'list' ){
        print $mf "$sname\t\t\t$sname\n";
    }else{die "choosing from mf|list for write mf!\n";}
    }
    close GP;
}

#=============
sub write_main{
#=============
    my ($handel,$STEP,$dir,$shell,$qopt ) = @_;
    my $middle;my $i=0;
    foreach my $shell (@$shell){
        my $mid_qsub=$qsub;
        $mid_qsub .= " $$qopt[$i] " if $qopt;
        $i++;
        $middle .= $qopt ? "$mid_qsub --input $shell 1>$shell.log 2>$shell.err\n" : "nohup $shell \n";
    }
    my $end = "date +\"\%D \%T: Finish $STEP\"";
    print $handel "#$STEP\ndate +\"\%D \%T: Start  $STEP\"\ncd $dir\n$middle$end\n\n";
}

#=============
sub write_step1{
#=============
    #outdir, shell dir, input data_list, mf
    my($step1_outdir,$step1_shell,$list,$mf)=@_;
    (-s $step1_shell) || `mkdir -p $step1_shell`;
    (-s $step1_outdir) || `mkdir -p $step1_outdir`;
    open(STEP1,">$step1_shell/step1.combine.qc.sh");
    open(LIST,">$step1_outdir/cleandata.list");
    open(OR,$list);
    while (<OR>) {
        chomp;
        my($sample,$fqs)=(split/\s+/)[0,1];
        (-d "$step1_outdir/$sample") || `mkdir -p $step1_outdir/$sample`;
        if($opt{type} eq 'pe'){
            my@fqs=split/,/,$fqs;
            print STEP1 
"cd $step1_outdir/$sample
mkdir -p 00.rawdata\n";
            if($opt{primer}){
                print STEP1 "#$primer_cut $fqs[0] $fqs[1] $opt{primer} $sample ./00.rawdata\n";       
		print STEP1 
"cutadapt  -b $opt{adapt}  -o 00.rawdata/$sample.raw_1.fq -p 00.rawdata/$sample.raw_2.fq $fqs[0] $fqs[1] >00.rawdata/cutadapt.log
gzip 00.rawdata/$sample.raw_1.fq -c >00.rawdata/$sample.raw_1.fq.gz
gzip 00.rawdata/$sample.raw_2.fq -c >00.rawdata/$sample.raw_2.fq.gz\n";
            }else{
                $fqs[0]=~/\.gz$/ ?
                print STEP1
"ln -fs $fqs[0] 00.rawdata/${sample}.raw_1.fq.gz
ln -fs $fqs[1] 00.rawdata/${sample}.raw_2.fq.gz\n":
                $fqs[0]=~/\.(fq|fastq)$/ ?
                print STEP1 
"gzip $fqs[0] -c > 00.rawdata/${sample}.raw_1.fq.gz
gzip $fqs[1] -c > 00.rawdata/${sample}.raw_2.fq.gz\n" :
                die "the raw data is not fq|fastq|fq.gz?\n";
            }

        print STEP1
"#$fqcheck -r 00.rawdata/${sample}.raw_1.fq.gz -c 00.rawdata/${sample}.raw_1.check & 
#$fqcheck -r 00.rawdata/${sample}.raw_2.fq.gz -c 00.rawdata/${sample}.raw_2.check &
$fqstat 00.rawdata/${sample}.raw_1.fq.gz,00.rawdata/${sample}.raw_2.fq.gz ./00.rawdata/ $sample.raw.fqcheck
#wait\n";
        $opt{cut} ? 
        print STEP1
"$read_cut 00.rawdata/${sample}.raw_1.fq.gz $r1_cut 00.rawdata/${sample}.cut$r1_cut.raw_1.fq.gz
$read_cut 00.rawdata/${sample}.raw_2.fq.gz $r2_cut 00.rawdata/${sample}.cut$r2_cut.raw_2.fq.gz
mkdir -p 01.combine
$flash 00.rawdata/${sample}.cut$r1_cut.raw_1.fq.gz 00.rawdata/${sample}.cut$r2_cut.raw_2.fq.gz -o 01.combine/$sample > 01.combine/$sample.flash.log
$fqsort 01.combine/$sample.extendedFrags.fastq  01.combine/$sample.extendedFrags.sort.fastq 300
mkdir -p 02.qc 
$split_libraries_fastq --sample_id $sample -o ./02.qc -i 01.combine/$sample.extendedFrags.sort.fastq -m $mf
":
        print STEP1
"mkdir -p 01.combine
$flash 00.rawdata/${sample}.raw_1.fq.gz 00.rawdata/${sample}.raw_2.fq.gz -o 01.combine/$sample > 01.combine/$sample.flash.log
$fqsort 01.combine/$sample.extendedFrags.fastq  01.combine/$sample.extendedFrags.sort.fastq 300
mkdir -p 02.qc 
$split_libraries_fastq --sample_id $sample -o ./02.qc -i 01.combine/$sample.extendedFrags.sort.fastq -m $mf
";
        }else{
            print STEP1 
"cd $step1_outdir/$sample
mkdir -p 00.rawdata\n";
        $fqs=~/\.gz$/ ?
        print STEP1 "ln -fs $fqs 00.rawdata/${sample}.raw.fq.gz\n" :
        $fqs=~/\.(fq|fastq)$/ ?
        print STEP1 "gzip $fqs -c > 00.rawdata/${sample}.raw.fq.gz\n" :
        die "the raw data is not fq|fastq|fq.gz?\n";
        print STEP1
"$fqcheck -r 00.rawdata/${sample}.raw.fq.gz -c 00.rawdata/${sample}.raw.check
mkdir -p 02.qc 
$split_libraries_fastq --sample_id $sample -o ./02.qc -i 00.rawdata/${sample}.raw.fq.gz -m $mf --phred_offset 33
";
        }
        if ($opt{region} =~/^16S/i) {
            print STEP1 
"mkdir -p 03.uchime
$usearch --db $gold -uchime_ref 02.qc/seqs.fna -strand plus -nonchimeras 03.uchime/seqs.nonchimeras.fna
perl $step1_lib/get.des.fq.pl 02.qc/seqs.fastq 03.uchime/seqs.nonchimeras.fna 03.uchime/seqs.nonchimeras.fastq
ln -fs 03.uchime/seqs.nonchimeras.fastq $sample.fastq
ln -fs 03.uchime/seqs.nonchimeras.fna $sample.fna
perl $step1_lib/len.draw.pl $sample.fna ./ $sample.tags.len
$ng_QC -i $sample.fastq -o 04.ng_QC
$fqstat $sample.fastq ./ $sample.fqcheck\n\n";
        }else{
            print STEP1 
"ln -fs 02.qc/seqs.fna $sample.fna
ln -fs 02.qc/seqs.fastq $sample.fastq
perl $step1_lib/len.draw.pl $sample.fna ./ $sample.tags.len
$ng_QC -i $sample.fastq -o 04.ng_QC
$fqstat $sample.fastq ./ $sample.fqcheck\n\n";
        }
        print LIST "$sample\t$step1_outdir/$sample/$sample.fna\n";
    }
    close LIST;
    close OR;
    close STEP1;
    open(STEP1,">$step1_shell/step2.qcstat.sh");
    print STEP1
"perl $step1_lib/qc_stat.pl $step1_outdir $mf  $step1_outdir/QCstat.xls
cd $step1_outdir
perl $Bin/../lib/00.Commbin/get.qc.report.pl ./ $mf ./\n";
    close STEP1;
    return( "$step1_outdir/cleandata.list","step1.combine.qc","step2.qcstat");
}

#=============
sub write_step2{
#=============
    #outdir, shell dir, input data_list, mf
    my($step2_outdir,$step2_shell,$list,$mf)=@_;
    (-s $step2_shell) || `mkdir -p $step2_shell`;
    (-s $step2_outdir) || `mkdir -p $step2_outdir`;
    open(OR,$list);
    my$total_fas;
    while (<OR>) {
        chomp;
        my@or=split/\s+/;
        $total_fas .=" $or[1]";
    }
    close OR;
    
    open(STEP1,">$step2_shell/step1.make.otu.sh");
    print STEP1
"#make otu table
mkdir -p $step2_outdir/01.Make_OTU
cd $step2_outdir/01.Make_OTU
cat  $total_fas > all.fna\n";
    if ($opt{otu} eq "uparse"){
        print STEP1
##$usearch11 -fastx_uniques all.fna -fastaout uniques.fasta -sizeout -minuniquesize 2
##$usearch11 -cluster_otus uniques.fasta -otus all_rep_set.fasta  -uparseout out.up -relabel OTU 
##$usearch11 -usearch_global all.fna -db all_rep_set.fasta -strand plus -id 0.97 -sample_delim _ -otutabout all_otu.txt
"#uparse64 
$uparse64 -derep_fulllength all.fna -output derep.fa -sizeout
$uparse64 -sortbysize derep.fa -output derep2.fa -minsize 2
$uparse64 -cluster_otus derep2.fa -otus outs1.fa
python $step2_lib/uparse_sc/fasta_number.py  outs1.fa OTU_ > all_rep_set.fasta
$uparse64 -usearch_global all.fna -db all_rep_set.fasta -strand plus -id 0.97 -uc map.uc
python $step2_lib/uparse_sc/uc2otutab.py map.uc > all_otu.txt
$assign_taxonomy -i all_rep_set.fasta -o rdp_assigned_taxonomy --c 0.8 --rdp_max_memory 18000
python $step2_lib/make_otu_table.py  all_otu.txt  rdp_assigned_taxonomy/all_rep_set_tax_assignments.txt  01.otu_table.txt
sed -i '1s/^/# Constructed from biom file\\n/' 01.otu_table.txt
perl $step2_lib/tags.stat.pl 01.otu_table.txt $mf 02.otu_table.new.txt Tags_stat.xls
biom convert -i 02.otu_table.new.txt -o 02.otu_table.new.biom --table-type \"OTU table\" --process-obs-metadata taxonomy --to-hdf5\n";
    }elsif($opt{otu} eq "usearch"){
        print STEP1
"$qiime/pick_otus.py -i all.fna -o ./01.pick.otus --denovo_otu_id_prefix OTU --threads 3
$qiime/pick_rep_set.py -i 01.pick.otus/all_otus.txt -f all.fna -o all_rep_set.fasta
$assign_taxonomy -i all_rep_set.fasta -o rdp_assigned_taxonomy --c 0.8 --rdp_max_memory 18000
$qiime/make_otu_table.py -i 01.pick.otus/all_otus.txt -t rdp_assigned_taxonomy/all_rep_set_tax_assignments.txt -o 01.otu_table.biom
biom convert -i 01.otu_table.biom -o 01.otu_table.txt --header-key taxonomy --to-tsv
perl $step2_lib/tags.stat.pl 01.otu_table.txt $mf 02.otu_table.new.txt Tags_stat.xls
biom convert -i 02.otu_table.new.txt -o 02.otu_table.new.biom --table-type \"OTU table\" --process-obs-metadata taxonomy --to-hdf5\n\n";
    }else{
        die "Otu method is not uparse and usearch\n\n";
    }
 
    print STEP1
"#get trees for otus
$qiime/align_seqs.py -i all_rep_set.fasta -t $core_set -o pynast_aligned_defaults/
$qiime/filter_alignment.py -i pynast_aligned_defaults/all_rep_set_aligned.fasta -o pynast_aligned_defaults/  
$qiime/make_phylogeny.py -i pynast_aligned_defaults/all_rep_set_aligned_pfiltered.fasta -o rep_set.tre\n\n" if $opt{region} =~ /^16s/i;

    print STEP1
"$qiime/sort_otu_table.py -i 02.otu_table.new.biom -o 03.sorted_otu_table.biom  -l $mf
biom summarize-table -i 03.sorted_otu_table.biom > 03.summarize.sorted.otu.table.txt
perl $step2_lib/rarefaction.pl 03.sorted_otu_table.biom 03.summarize.sorted.otu.table.txt 04.otu_table.even.biom min single
biom convert -i 04.otu_table.even.biom -o 04.otu_table_with_taxonomy.txt --header-key taxonomy --to-tsv
mkdir -p $step2_outdir/02.Make_Tables
cd $step2_outdir/02.Make_Tables
$qiime/summarize_taxa.py -i $step2_outdir/01.Make_OTU/04.otu_table.even.biom -o Absolute --suppress_biom_table_output -a 
$qiime/summarize_taxa.py -i $step2_outdir/01.Make_OTU/04.otu_table.even.biom -o Rel --suppress_biom_table_output
perl $step2_lib/trans.relative.pl Rel Relative tax_table\n\n";
    if ($opt{group}) {
        print STEP1
"#make tables for group
$avg_group $step2_outdir/02.Make_Tables/Relative $mf $step2_outdir/02.Make_Tables/Rel.group
$sort_max $step2_outdir/02.Make_Tables/Rel.group $step2_outdir/02.Make_Tables/Relative.group\n";
    }
    close STEP1;

    open(STEP3,">$step2_shell/step2.make.stats.sh");
    if( $opt{mf}){
        print STEP3 
"#top & heatmap for samples
perl $step4_lib/group.list.pl $mf $step2_outdir/group.list
$Rscript $top_heatmap --indir $step2_outdir/02.Make_Tables/Relative --outdir $step2_outdir/03.Tax_Stats --prefix tax_table  --top $opt{top} --heatmap $opt{heatmap} --group $step2_outdir/group.list\n\n" ;
    }else{
        print STEP3 
"#top & heatmap for samples
$Rscript $top_heatmap --indir $step2_outdir/02.Make_Tables/Relative --outdir $step2_outdir/03.Tax_Stats --prefix tax_table  --top $opt{top} --heatmap $opt{heatmap}\n\n";
    }

    print STEP3 
"#top & heatmap for groups
$Rscript $step2_lib/top_heatmap.r  --indir $step2_outdir/02.Make_Tables/Relative.group --outdir $step2_outdir/03.Tax_Stats.group --prefix tax_table  --top $opt{top} --heatmap $opt{heatmap}\n\n" if $opt{group} && $opt{group} =~/1/;

    print STEP3
"#Kurskal-Wallis for taxonomy
mkdir -p $step2_outdir/03.Tax_Stats.group/Kurskal_Wallis/01.testResult
cd $step2_outdir/03.Tax_Stats.group/Kurskal_Wallis/01.testResult
$kur_wal $step2_outdir/02.Make_Tables/Relative/tax_table.p.txt $mf phylum.KW 0
$kur_wal $step2_outdir/02.Make_Tables/Relative/tax_table.c.txt $mf class.KW 0
$kur_wal $step2_outdir/02.Make_Tables/Relative/tax_table.o.txt $mf order.KW 0
$kur_wal $step2_outdir/02.Make_Tables/Relative/tax_table.f.txt $mf family.KW 0
$kur_wal $step2_outdir/02.Make_Tables/Relative/tax_table.g.txt $mf genus.KW 0
$Rscript $step2_lib/top_heatmap.KW.r --indir $step2_outdir/03.Tax_Stats.group/Kurskal_Wallis/01.testResult --outdir $step2_outdir/03.Tax_Stats.group/Kurskal_Wallis/02.top.Heatmap --top $opt{top} --heatmap $opt{heatmap} --group $step2_outdir/group.list\n\n" if $opt{group} && $opt{group}=~/2/;

    print STEP3 
"#LefSe for groups
mkdir -p $step2_outdir/03.Tax_Stats.group/LefSe
cd $step2_outdir/03.Tax_Stats.group/LefSe
$otu2lefse $step2_outdir/02.Make_Tables/Relative/ $mf LefSe.input.txt
$lefse/format_input.py LefSe.input.txt LefSe.input.in -c 1 -o 1000000
$lefse/run_lefse.py LefSe.input.in LefSe.output -l $opt{lda}
$lefse/plot_res.py --left_space 0.3 LefSe.output LefSe.barplot.pdf --format pdf
/usr/bin/convert -density 300 LefSe.barplot.pdf LefSe.barplot.png
$lefse/plot_cladogram.py LefSe.output LefSe.cladogram.pdf --format pdf
/usr/bin/convert -density 300 LefSe.cladogram.pdf LefSe.cladogram.png
mkdir ./plot_features/
$lefse/plot_features.py LefSe.input.in LefSe.output  ./plot_features/ --format pdf\n\n" if $opt{group} && $opt{group}=~/3/; 
    
    if ($opt{venn}) {
        -s $opt{venn} || print "Warnings:venn list is not exists\n";
        $opt{venn}=abs_path($opt{venn});
        print STEP3 
"#venn for samples
mkdir -p $step2_outdir/03.Tax_Stats/Venn
cd $step2_outdir/03.Tax_Stats/Venn
perl $step2_lib/format.otu.table.pl $step2_outdir/01.Make_OTU/04.otu_table_with_taxonomy.txt otu_table.txt
perl $step2_lib/format.table.pl $opt{venn} venn.list
$Rscript  $step2_lib/venn.R --input $step2_outdir/03.Tax_Stats/Venn/otu_table.txt --venn $step2_outdir/03.Tax_Stats/Venn/venn.list\n\n";
    }

    if($opt{venn_group}){
        -s $opt{venn_group} || print "Warnings:venn group list is not exists\n";
        $opt{venn_group}=abs_path($opt{venn_group});
        print STEP3
"#venn for groups
mkdir -p $step2_outdir/03.Tax_Stats/venn_group
cd $step2_outdir/03.Tax_Stats/venn_group
perl $step2_lib/format.table.pl $opt{venn_group} venn.group.list
$qiime/collapse_samples.py -m $mf -b $step2_outdir/01.Make_OTU/04.otu_table.even.biom --output_biom_fp $step2_outdir/03.Tax_Stats/venn_group/Description_otu_table.biom --output_mapping_fp $step2_outdir/03.Tax_Stats/venn_group/Description_map.txt --collapse_fields 'Description'
biom convert -i $step2_outdir/03.Tax_Stats/venn_group/Description_otu_table.biom -o $step2_outdir/03.Tax_Stats/venn_group/Description_otu_table.txt --header-key taxonomy --to-tsv
perl $step2_lib/format.otu.table.pl $step2_outdir/03.Tax_Stats/venn_group/Description_otu_table.txt otu_table.txt
$Rscript  $step2_lib/venn.R --input $step2_outdir/03.Tax_Stats/venn_group/otu_table.txt --venn $step2_outdir/03.Tax_Stats/venn_group/venn.group.list\n\n";
    }

    close STEP3;
    return("$step2_outdir/01.Make_OTU/04.otu_table_with_taxonomy.txt","$step2_outdir/01.Make_OTU/04.otu_table.even.biom","$step2_outdir/01.Make_OTU/rep_set.tre","$step2_outdir/01.Make_OTU/03.summarize.sorted.otu.table.txt","step1.make.otu","step2.make.stats");
}
#=============
sub write_step3{
#=============
    my($otu_table,$otu_biom,$tree,$summarize,$mf,$step3_shell,$step3_outdir)=@_;
    (-s $step3_shell) || `mkdir -p $step3_shell`;
    (-s $step3_outdir) || `mkdir -p $step3_outdir`;
    open(STEP,">$step3_shell/alpha.diversity.sh");
    print STEP
"cd $step3_outdir
perl $step3_lib/group.list.pl $mf  $step3_outdir/alpha.mf

#$qiime/alpha_rarefaction.py -i $otu_biom -m $step3_outdir/alpha.mf -p $step3_lib/beta_params.txt -t $tree -o 01.alpha_diversity -n 6
#Alpha rarefaction command
perl $step2_lib/rarefaction.pl $otu_biom $summarize ./01.rarefaction/ min multiple
# Alpha diversity on rarefied OTU tables command
$qiime/alpha_diversity.py -i ./01.rarefaction/ -o ./02.alpha_div/  -t $tree -m observed_otus,shannon,chao1
# Collate alpha command
$qiime/collate_alpha.py -i ./02.alpha_div/ -o ./03.alpha_div_collated/
# Removing intermediate files command
#rm -r ./01.rarefaction/ ./02.alpha_div/
# Rarefaction plot: All metrics command
$qiime/make_rarefaction_plots.py -i ./03.alpha_div_collated/ -m $step3_outdir/alpha.mf -o ./04.alpha_rarefaction_plots/
mkdir -p 05.alpha_plots/
perl $step3_lib/alpha.input.pl ./03.alpha_div_collated/chao1.txt $step3_outdir/05.alpha_plots/chao1.avg.txt
$Rscript $rare --infile $step3_outdir/05.alpha_plots/chao1.avg.txt --outdir  $step3_outdir/05.alpha_plots --out Chao1_index
perl $step3_lib/alpha.input.pl ./03.alpha_div_collated/shannon.txt $step3_outdir/05.alpha_plots/shannon.avg.txt
$Rscript $rare  --infile $step3_outdir/05.alpha_plots/shannon.avg.txt --outdir  $step3_outdir/05.alpha_plots --out Shannon_index
perl $step3_lib/alpha.input.pl ./03.alpha_div_collated/observed_otus.txt $step3_outdir/05.alpha_plots/observed_otus.avg.txt
$Rscript $rare  --infile $step3_outdir/05.alpha_plots/observed_otus.avg.txt --outdir  $step3_outdir/05.alpha_plots --out Observed_otus

perl $step3_lib/format.otu.table.pl $otu_table otu_table.txt
$Rscript $rank --infile otu_table.txt --out rank_abundance.pdf
convert -density 300 rank_abundance.pdf rank_abundance.png
$qiime/alpha_diversity.py -i $otu_biom -m observed_otus,chao1,shannon,simpson,goods_coverage -t $tree -o alpha_diversity.txt
perl $step3_lib/alpha_index_format.pl alpha_diversity.txt alpha_diversity_index.xls
rm alpha_diversity.txt\n\n";

    if ($opt{group} && $opt{group} =~/4/) {
        print STEP 
"#Kurskal-Wallis for alpha_diversity
perl $step3_lib/trans.table.pl alpha_diversity_index.xls > alpha_diversity_index.xls.trans
$kw_alpha alpha_diversity_index.xls.trans $mf $step3_outdir/06.KW.ByGroups alpha.KW
convert -density 300 06.KW.ByGroups/observed_otus.violin.pdf 06.KW.ByGroups/observed_otus.violin.png
convert -density 300 06.KW.ByGroups/chao1.violin.pdf 06.KW.ByGroups/chao1.violin.png
convert -density 300 06.KW.ByGroups/shannon.violin.pdf 06.KW.ByGroups/shannon.violin.png
convert -density 300 06.KW.ByGroups/simpson.violin.pdf 06.KW.ByGroups/simpson.violin.png\n";
    }

    close STEP;
    return("alpha.diversity");
}
#=============
sub write_step4{
#=============
    my($otu_biom,$tree,$mf,$outdir,$shell)=@_;
    (-s $outdir) || `mkdir -p $outdir`;
    (-s $shell) || `mkdir -p $shell`;
    open(STEP,">$shell/beta_diversity.sh");
    print STEP
"cd $outdir
# Beta Diversity (weighted_unifrac) command
$qiime/beta_diversity.py -i $otu_biom -o ./01.1.bdiv_even/ --metrics weighted_unifrac  -t $tree

# Rename distance matrix (weighted_unifrac) command
mv ./01.1.bdiv_even/weighted_unifrac_04.otu_table.even.txt ./01.1.bdiv_even/weighted_unifrac_dm.txt

# Principal coordinates (weighted_unifrac) command
$qiime/principal_coordinates.py -i ./01.1.bdiv_even/weighted_unifrac_dm.txt -o ./01.1.bdiv_even/weighted_unifrac_pc.txt

# Make emperor plots, (weighted_unifrac) command
$qiime/make_emperor.py -i ./01.1.bdiv_even/weighted_unifrac_pc.txt -o 01.2.weighted_unifrac_emperor_pcoa_plot/ -m $mf

# Make pcoa 2d plot, (weighted_unifrac) command
perl $step4_lib/group.list.pl $mf  $outdir/group.list
$Rscript $pcoa --indm $outdir/01.1.bdiv_even/weighted_unifrac_dm.txt --outdir $outdir//01.3.weighted_2d_plot --group $outdir/group.list
convert -density 300 $outdir/01.3.weighted_2d_plot/Dendrogram.pdf $outdir/01.3.weighted_2d_plot/Dendrogram.png
convert -density 300 $outdir/01.3.weighted_2d_plot/PCoA12-2.pdf $outdir/01.3.weighted_2d_plot/PCoA12-2.png 
convert -density 300 $outdir/01.3.weighted_2d_plot/PCoA12.pdf $outdir/01.3.weighted_2d_plot/PCoA12.png 

# Beta Diversity (unweighted_unifrac) command
$qiime/beta_diversity.py -i $otu_biom -o ./02.1.bdiv_even/ --metrics unweighted_unifrac  -t $tree

# Rename distance matrix (unweighted_unifrac) command
mv ./02.1.bdiv_even/unweighted_unifrac_04.otu_table.even.txt ./02.1.bdiv_even/unweighted_unifrac_dm.txt

# Principal coordinates (unweighted_unifrac) command
$qiime/principal_coordinates.py -i ./02.1.bdiv_even/unweighted_unifrac_dm.txt -o ./02.1.bdiv_even/unweighted_unifrac_pc.txt

# Make emperor plots, unweighted_unifrac) command
$qiime/make_emperor.py -i ./02.1.bdiv_even/unweighted_unifrac_pc.txt -o ./02.2.unweighted_unifrac_emperor_pcoa_plot/ -m $mf

# Make pcoa 2d plot, (unweighted_unifrac) command
$Rscript $pcoa --indm $outdir/02.1.bdiv_even/unweighted_unifrac_dm.txt --outdir $outdir/02.3.unweighted_2d_plot --group $outdir/group.list
convert -density 300 $outdir/02.3.unweighted_2d_plot/Dendrogram.pdf $outdir/02.3.unweighted_2d_plot/Dendrogram.png
convert -density 300 $outdir/02.3.unweighted_2d_plot/PCoA12-2.pdf $outdir/02.3.unweighted_2d_plot/PCoA12-2.png
convert -density 300 $outdir/02.3.unweighted_2d_plot/PCoA12.pdf $outdir/02.3.unweighted_2d_plot/PCoA12.png\n";
    close STEP;
    return"beta_diversity";
}
#====================================================================================================================
