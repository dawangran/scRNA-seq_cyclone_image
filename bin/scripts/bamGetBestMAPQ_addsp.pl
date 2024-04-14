#!/usr/bin/perl -w
use strict;

my $header=shift;
my $file=shift;
my $mapq_the=shift;


open IN,"$header";
while(<IN>){
    print $_;
}close IN;

my %hash;
open IN,"samtools view $file|";
# open IN,"$file";
while(<IN>){
    chomp;
    my @tmp=split /\t/;
    my $l1=$tmp[0];
    my $mapq=$tmp[4];

    shift @tmp;
    my $line=join("\t",@tmp);

    my $readid=(split /#/,$l1)[-1];
    my $cbumi=(split /#/,$l1)[0];

    my $umi=(split /_/,$cbumi)[-1];
    my $barcode=(split /_/,$cbumi)[0];

    my $line_new=$readid."\t".$line."\tUR:Z:".$umi."\tCB:Z:".$barcode;

    if($mapq >= $mapq_the && $hash{$readid}){
        my $mapqnow=(split /\t/,$hash{$readid})[4];
        if($mapq > $mapqnow){
            $hash{$readid}=$line_new;
        }
    }elsif($mapq >= $mapq_the ){
        $hash{$readid}=$line_new;
    }
}close IN;

for my $key (keys %hash){
    
    print $hash{$key};
    my $sp=(split/\t/,$hash{$key})[2];
    if($sp =~ /GRCh38_/){
        print "\tSP:Z:Human\n";
    }elsif($sp =~ /mm10_/){
        print "\tSP:Z:Mouse\n";
    }else{
        print "\n";
    }

}
