#!/usr/bin/perl
use strict;
use Getopt::Long;

my $maf;
my $maxmaf;
GetOptions(
        "minaf:s" => \$maf,
	"maxaf:s" => \$maxmaf
);

$maf ||= 0.0001;
$maxmaf ||=0.5;
#CHR POS REF ALT ChinaMAP gnomAD-EAS EAS GT
my $num=0;
open(OT1,">tmp.af");
open(OT2,">tmp.tped");
open(OT3,">mogan.s.bed");
while(<>){
         chomp;
         my @t=(split/\s+/,$_,7);
	 my $ref=$t[2];
	 my $alt=$t[3];
         my $af=$t[4];
         $af=$t[5] if($af=='.');
	 $af=$maf if($af=='.');
	 next if($af < $maf || $af > $maxmaf);
	 $num++;
	 my $af1=1-$af;
         #$af=0.001 if($af=='.' || $af==0 || $af<0.001);
	 #$af=0.999 if($af>0.999);
	 $t[-1]=~s/\/|\|/ /g;
	 $t[-1]=~s/1/2/g;
	 $t[-1]=~s/0/1/g;
	 $t[-1]=~s/\./0/g;
	 print OT1 "set markers $num allele freq $af1 $af\n";
	 print OT2 "$t[0] maker$num 0 $t[1] $t[-1]\n";
	 print OT3 "$t[0]\t",$t[1]-1,"\t$t[1]\n";
}
print OT1 "\nset marker data $num\n";
