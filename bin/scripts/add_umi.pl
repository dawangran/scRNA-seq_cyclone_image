#!/usr/bin/perl -w
use strict;

my $fq=shift;


open IN,"$fq";
while(my $line1=<IN>){
    chomp$line1;
    my $line2=<IN>;chomp$line2;
    my $line3=<IN>;chomp$line3;
    my $line4=<IN>;chomp$line4;
    
    if(length($line2) > 20){
        my $umi=substr($line2,0,10);
        my $seq=substr($line2,10);
        my $qual=substr($line4,10);

        print $line1."|||UR:Z:".$umi."\n";
        print $seq."\n";
        print "+\n";
        print $qual."\n";
    }
   



}close IN;
