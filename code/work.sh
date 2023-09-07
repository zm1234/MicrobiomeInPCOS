#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N main.16s
#$ -q NGS
#1) run Data Control
date +"%D %T: Start  1) run Data Control"
cd Phase1/01.standeranalysis/Shell/detail/01.Combine_QC/
perl bin/../lib/00.Commbin/super.qsub.pl   --prefix combine_qc --dir step1.combine.qc   --input step1.combine.qc.sh 1>step1.combine.qc.sh.log 2>step1.combine.qc.sh.err
perl bin/../lib/00.Commbin/super.qsub.pl   --prefix qc_stat --dir step2.qcstat   --input step2.qcstat.sh 1>step2.qcstat.sh.log 2>step2.qcstat.sh.err
date +"%D %T: Finish 1) run Data Control"

#2) run OTU analysis
date +"%D %T: Start  2) run OTU analysis"
cd Phase1/01.standeranalysis/Shell/detail/02.OTU_analysis/
perl bin/../lib/00.Commbin/super.qsub.pl   --prefix make_otu --dir step1.make.otu --parallel 1 --ppn 12   --input step1.make.otu.sh 1>step1.make.otu.sh.log 2>step1.make.otu.sh.err
perl bin/../lib/00.Commbin/super.qsub.pl   --prefix otu_stat --dir step2.make.stats   --input step2.make.stats.sh 1>step2.make.stats.sh.log 2>step2.make.stats.sh.err
date +"%D %T: Finish 2) run OTU analysis"

#3) run alpha diveristy
date +"%D %T: Start  3) run alpha diveristy"
cd Phase1/01.standeranalysis/Shell/detail/03.Alpha_diversity/
perl bin/../lib/00.Commbin/super.qsub.pl   --prefix alpha --dir alpha.diversity --parallel 1   --input alpha.diversity.sh 1>alpha.diversity.sh.log 2>alpha.diversity.sh.err
date +"%D %T: Finish 3) run alpha diveristy"

#4) run beta diveristy
date +"%D %T: Start  4) run beta diveristy"
cd Phase1/01.standeranalysis/Shell/detail/04.Beta_diveristy/
perl bin/../lib/00.Commbin/super.qsub.pl   --prefix beta --dir beta_diversity --parallel 1   --input beta_diversity.sh 1>beta_diversity.sh.log 2>beta_diversity.sh.err
date +"%D %T: Finish 4) run beta diveristy"

#5) get result&report
date +"%D %T: Start  5) get result&report"
cd Phase1/01.standeranalysis/Shell/detail
nohup perl bin/../lib/00.Commbin/get.result.pl  Phase1/01.standeranalysis Phase1/01.standeranalysis/Shell/samples.mf Phase1/01.standeranalysis/Shell/detail/05.result_report/ Phase1/01.standeranalysis/
perl bin/../lib/00.Commbin/get.report.pl  Phase1/01.standeranalysis/result Phase1/01.standeranalysis 
date +"%D %T: Finish 5) get result&report"

