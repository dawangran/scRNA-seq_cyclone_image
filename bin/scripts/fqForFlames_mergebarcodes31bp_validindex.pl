#!/usr/bin/perl -w
use strict;

my $barcodefile=shift;
my $infile=shift;
my $valid=shift; # 组合
my $num=shift; # 区分lib名字

$/=">";my %hash;
open IN,"$barcodefile";
<IN>;
while(<IN>){
    chomp;
    my @tmp=split /\n/;
    $hash{$tmp[0]}=$tmp[-1];
}close IN;

$/="\n";

my %hashv;
open IN,"$valid";
while(<IN>){
    chomp;
    my @tmp=split /\t/;
    my $a=substr($tmp[0],0,10);
    my $b=substr($tmp[0],10,10);
    # my $new=$a."CCTTCC".$b."CGATG";
    my $new=$a."CCTTCC".$b;
    $tmp[1] =~ s/_//g;
    $hashv{$new}=$tmp[1];
}close IN;






open IN,$infile;
while(my $line1=<IN>){
    chomp$line1;
    my $line2=<IN>;chomp$line2;
    my $line3=<IN>;chomp$line3;
    my $line4=<IN>;chomp$line4;

    my @tmp=split /\|\|\|/, $line1;
    $tmp[0] =~ s/\@/\#/;
    $tmp[1]=(split /:/,$tmp[1])[-1];
    $tmp[2]=(split /:/,$tmp[2])[-1];
    my $barcode=$hash{$tmp[1]};
    if($hashv{$barcode}){
        print "@".$hashv{$barcode}."_".$tmp[2].$tmp[0]."lib".$num."\n";
        print $line2."\n".$line3."\n".$line4."\n";
    }

    



}close IN;
