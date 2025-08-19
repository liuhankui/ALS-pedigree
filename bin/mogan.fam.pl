#!/usr/bin/perl
use strict;

my $num=0;
my $o="";
while(<>){
         chomp;
         my @t=split;
         $num++;
         $o .="$t[1] $t[2] $t[3] $t[4] $t[5]\n";
}

print "input pedigree size $num\n";
print "input pedigree record names 3 integers 2\n";
print "input pedigree record trait 1 integers 2\n";
print "*****\n";
print "$o";
