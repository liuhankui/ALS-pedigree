#!/usr/bin/perl
use strict;
use Getopt::Long;

my ($fam);
GetOptions(
        "f:s" => \$fam,
);

my %error=(
         "0/0 0/1 1/1"=>1,
         "0/0 1/1 0/1"=>1,
         "0/0 1/1 1/1"=>1,
         "0/0 0/0 1/1"=>1,
         "0/0 1/1 0/0"=>1,

         "0/0 ./. 1/1"=>1,
         "0/0 1/1 ./."=>1,

         "0/1 0/0 0/0"=>1,
         "0/1 1/1 1/1"=>1,

         "1/1 0/0 0/0"=>1,
         "1/1 0/1 0/0"=>1,
         "1/1 0/0 0/1"=>1,
         "1/1 0/0 1/1"=>1,
         "1/1 1/1 0/0"=>1,

         "1/1 0/0 ./."=>1,
         "1/1 ./. 0/0"=>1,
);

my %fam;
open(FAM,"$fam");
while(<FAM>){
         chomp;
         my @t=split;
         next if($t[2] eq '0' || $t[3] eq '0');
         $fam{$t[1]}[0]=$t[2];
         $fam{$t[1]}[1]=$t[3];
}

my @id;
my %order;
open(EO,">snv.md.error");
while(<>){
         chomp;
         if(/#/){
                 if(/#CH/){
                         @id=(split/\t/,(split/\t/,$_,10)[-1]);
                         for(my $i=0;$i<=$#id;$i++){
                                 $order{$i}=$id[$i];
                         }
                 }
		 print "$_\n";
                 next;
         }

         my @t=(split/\t/,$_,10);
         $t[-1]=~s/\|/\//g;
         my @g=(split/\t/,$t[-1]);
         my %geno;

         for(my $i=0;$i<=$#id;$i++){
                 my $id=$id[$i];
                 $geno{$id}=(split/\:/,$g[$i])[0];
         }

         my $err=0;
         for(my $i=0;$i<=$#id;$i++){
                 my $id=$id[$i];
                 next if(!exists($fam{$id}));
                 my $fid=$fam{$id}[0];
                 my $mid=$fam{$id}[1];
		 my $cgeno=$geno{$id} || './.';
                 my $fgeno=$geno{$fid} || './.';
                 my $mgeno=$geno{$mid} || './.';
                 my $fam_geno="$cgeno $fgeno $mgeno";
                 if(exists($error{$fam_geno})){
		  $err=1;
		  print EO "error: $fam_geno $id $fid $mid $t[1] \n";
		 }else{
#                  print EO "corrt: $fam_geno $id $fid $mid $t[1]\n";
		 }
		 
         }
         next if($err);
         print "$_\n";
}
close EO;
