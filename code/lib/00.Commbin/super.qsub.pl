use strict;
use warnings;

use Cwd qw(abs_path);
use FindBin qw($Bin);
use Getopt::Long;
use POSIX; 

#set default paramers
my %opt = (
	"split",'\n\n',"dir","./Shell","prefix","work","ppn",1,"nodes",1,"query","NGS","time",120,
);

#get options
GetOptions(
    \%opt,
    "input:s","split:s","dir:s","prefix:s","ppn:s","nodes:s","query:s","parallel:s","time:s","notrun",
    #"locate",
);

$opt{input} && -s $opt{input} || die
"Usage: perl $0 --input *.sh 

	--input  	[file]   input shell file for qsub
	--split  	[str]    split for qsub, default=\\n\\n
	--dir 	 	[str]    work directory for shells, default=./Shell
	--prefix 	[str]    shell prefix for split, default=work
	--ppn 		[str]    ppn' number for qsub, default=1 
	--nodes		[str]	 nodes' number for qsub, default=1
	--query		[str]	 query for qsub, default=NGS
	--parallel      [str]    the most shell number to parallel, default not set 
	--time          [str]    time for checking stats. default=120(s)
	--locate	         run the shell by locate, default not set（！not used yet）
	--notrun  		 just get pbs, not qsub, default not set 

";
###################################################################################################
#set options
-s $opt{dir} || `mkdir -p $opt{dir}`;
$opt{dir}=abs_path($opt{dir});
(-s "$opt{dir}/Shells$$") || `mkdir -p $opt{dir}/Shells$$`;

###################################################################################################
#main scripts
my $shells=`cat $opt{input}`;
chomp$shells;
my @shells=split/$opt{split}/,$shells;

#set for split num
my $split_num;
if($opt{parallel}){
	$split_num=int(@shells/$opt{parallel});
	$split_num++ if @shells%$opt{parallel};
}else{
	$split_num=1;
	$opt{parallel}=@shells;
}

my $i=1;my $m=1; #i for shell number, m for write number
open(PBS,">$opt{dir}/Shells$$/$opt{prefix}$m.pbs");
##PBS -l nodes=$opt{nodes}:ppn=$opt{ppn}
print PBS 
"#!/bin/bash
#\$ -cwd
#\$ -S /bin/bash
#\$ -N $opt{prefix}$m
#\$ -pe mpi $opt{ppn}
#\$ -q $opt{query}\n";
if (! $opt{locate}  && ! $opt{notrun}) {
	open(OUT,">$opt{dir}/qsub.pid.$$.out");
	open(QDEL,">$opt{dir}/qdel.$$.sh");
	print QDEL "qdel ";
}elsif($opt{locate}){
	print PBS "echo \"\$0 \$\$\">> $opt{dir}/locate.pid.out\n";
}

my@pids;
foreach my $sh (@shells){
	print PBS "$sh\n\n";
	if( $i%$split_num == 0){
		#the last writing
		print PBS "date|awk '{print \$0,\"finish $opt{prefix}$m.pbs\"}' >> $opt{dir}/finish.$$.out\n\n";
		close PBS;
		my $pid=&execu_shell($m);
		push @pids,$pid if $pid;
		#next writing
		$m++;
		last if $m > $opt{parallel} || $m > @shells;
		open(PBS,">$opt{dir}/Shells$$/$opt{prefix}$m.pbs");
##PBS -l nodes=$opt{nodes}:ppn=$opt{ppn}
		print PBS 
"#!/bin/bash
#\$ -cwd
#\$ -S /bin/bash
#\$ -N $opt{prefix}$m
#\$ -pe mpi $opt{ppn}
#\$ -q $opt{query}\n";
		print PBS "echo \"\$0 \$\$>> $opt{dir}/locate.pid.out\n" if $opt{locate};
	}
	$i++;
}

if ($opt{parallel} < @shells && $m < $opt{parallel}) {
	print PBS "date|awk '{print \$0,\"finish $opt{prefix}$m.pbs\"}' >> $opt{dir}/finish.$$.out\n\n";
	close PBS;
	my $pid=&execu_shell($m);
	push @pids,$pid if $pid;
}

close OUT if ! $opt{locate}  && ! $opt{notrun};
close QDEL if ! $opt{notrun};

if ($opt{time} && ! $opt{locate}  && ! $opt{notrun} ) {
	my $outf="$opt{dir}/qsub.stat.$$.out";
	&check_stat(\@pids,$outf);
}



###################################################################################################
###################################################################################################
##=====================================##
##            SUB FUNCTIONs            ##
##=====================================##
##############
#give a number, qsub or locate
sub execu_shell{
##############
	my($m,)=@_;
	my $prog="$opt{prefix}$m.pbs";
	my $pid;
		if (! $opt{locate}  && ! $opt{notrun}) {
			$pid=`cd $opt{dir}/Shells$$
				qsub -cwd $prog `;
			chomp$pid;
			$pid=$1 if $pid=~/Your job (\d+) .*has been submitted/;
			print OUT "$pid\t$prog\n";
			print QDEL "$pid ";
		}elsif($opt{locate}){
			sleep(3+rand(5));
			system"cd $opt{dir}/Shells$$;bg sh $prog & ";
		}
	return $pid;
}
##############
#give pids, checking for stat
sub check_stat{
##############
	my ($pids,$out)= @_;
	my%check;
	map{$check{$_}=1;} @$pids;
	chomp(my $user=`whoami`);
	open(OUTS,">$out");
	while (1) {
		my @stats=split/\n/,`qstat -u $user`;
		my $flag;
		foreach my $stat (@stats){
			$stat=~s/^\s+//;
			my@infos=split/\s+/,$stat;
			print OUTS `date`."$stat\n" if $infos[0] && $infos[0]=~/^(Job|--)/;
			print OUTS "$stat\n" if $infos[0] && $check{$infos[0]} ;
			$flag = 1 if $infos[0] && $check{$infos[0]};
		}
		last if !$flag;
		print OUTS "\n";
		sleep$opt{time};
	}
	close OUTS;
}
