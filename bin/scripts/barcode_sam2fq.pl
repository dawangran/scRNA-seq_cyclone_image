#!/usr/bin/perl -w
use strict;

my $sam=shift;
my $barcodefile=shift;


# $/=">";my %hash;
# open IN,"$barcodefile";
# <IN>;
# while(<IN>){
#     chomp;
#     my @tmp=split /\n/;
#     $hash{$tmp[0]}=$tmp[-1];
# }close IN;

my %hh;
$/="\n";
open IN,"$sam";
while(<IN>){
    chomp;
    my @tmp=split/\s+/;
    my $length=length($tmp[9]);
    if(!$hh{$tmp[0]}){
        my $count_m = $tmp[5] =~ tr/M/M/; 
        # print $tmp[5]."\n";
        if($count_m == 1 ){
            my $first=(split /M/,$tmp[5])[0];

            my $alpha_num=$first =~ tr/[A-Z]/[A-Z]/;
            if($alpha_num == 0){
                my $numM=$first;
                my $index_end=$numM-1;
                my $seq=substr($tmp[9],-($length-$index_end - 1));
                my $qual=substr($tmp[10],-($length-$index_end - 1));

                print "\@".$tmp[0]."|||CB:Z:".$tmp[2]."\n";
                print $seq."\n";
                print "+\n";
                print $qual."\n";

            }else{
                my $numS=(split /[A-Z]/,$first)[0];
                my $numM=(split /[A-Z]/,$first)[1];
                my $index_start=$numS;
                my $index_end=$numS+$numM-1;
                my $seq=substr($tmp[9],-($length-$index_end - 1));
                my $qual=substr($tmp[10],-($length-$index_end - 1));

                print "\@".$tmp[0]."|||CB:Z:".$tmp[2]."\n";
                print $seq."\n";
                print "+\n";
                print $qual."\n";
            }
        }else{
            my @tmp1=(split /M/,$tmp[5]);
            pop @tmp1;
            my $tline=join("N",@tmp1);
            my @num=split /[A-Z]/,$tline;

            my $total=0;
            map{$total+=$_} @num;

            
            # print $total;
            my $seq=substr($tmp[9],-($length-$total));
            my $qual=substr($tmp[10],-($length-$total));
            
            print "\@".$tmp[0]."|||CB:Z:".$tmp[2]."\n";
            print $seq."\n";
            print "+\n";
            print $qual."\n";

        }
        $hh{$tmp[0]}=1;
    }
    
    




}close IN;