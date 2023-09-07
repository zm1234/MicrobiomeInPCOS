#!/usr/bin/perl -w
use strict;
use PerlIO::gzip;
die"usage: perl <in.fq> <in.fa> <out.fa>\n" unless @ARGV==3;
my ($infq,$infa,$outfq) = @ARGV;

my %check;
open FA,$infa || die$!;
while (<FA>) {
    if(/>(\S+)/){
            $check{$1}=1;
    }
}
close FA;

$infq=~/\.gz$/ ? open(FQ,"<:gzip",$infq) : open(FQ,$infq);
$outfq=~/\.gz$/ ? open(OUT,">:gizp",$outfq) : open(OUT,">",$outfq);
while(my $head = <FQ>){
    my $id =$1 if $head=~/\@(\S+)/;
    if($check{$id}){
        $head .= (<FQ> . <FQ> . <FQ>);
        print OUT $head;
    }else{
        <FQ>;<FQ>;<FQ>;
    }
}
close FQ;
close OUT;