# ALS-pedigree

# WGS variant calling (SNV,INDEL,SV, CNV,STR) 

```
# Alignment
fastp -q 20 -u 10 -n 5 --in1 fq1.gz --in2 fq2.gz --out1 clean.fq1.gz --out2 clean.fq2.gz
bwa mem -t 4 -R '@RG\tID:id\tPL:illumina\tPU:id\tLB:sample\tSM:id\tCN:BGI' GRCh38.fa clean.fq1.gz clean.fq2.gz|samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 -d discordants.sam -s splitters.sam|samtools view -Sb -|samtools sort - -O CRAM -o sort.cram --reference GRCh38.fa
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" discordants.sam|sambamba view -S -f bam /dev/stdin|samtools sort - -O CRAM -o discordants.cram --reference GRCh38.fa
gawk '{ if ($0~"^@") { print; next } else { $10="*"; $11="*"; print } }' OFS="\t" splitters.sam|sambamba view -S -f bam -l 0 /dev/stdin|samtools sort - -O CRAM -o splitters.cram --reference GRCh38.fa
# QC
mosdepth qc sort.cram -f GRCh38.fa --fast-mode --no-per-base --by 1000000 --thresholds 1,2,4,10
# SNV/INDEL calling
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar HaplotypeCaller --QUIET true -R GRCh38.fa -I sort.cram -O gatk.gvcf.gz -ERC GVCF -A ClippingRankSumTest -A LikelihoodRankSumTest -A MappingQualityZero
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar GenotypeGVCFs --QUIET true -R GRCh38.fa --variant gatk.gvcf.gz -O gatk.vcf.gz
# SV calling
lumpyexpress -P -T ./tmp -R GRCh38.fa -B sort.cram -S splitters.bam -D discordants.bam -x exclude.bed -o sv.vcf && svtyper -B sort.cram -T GRCh38.fa -i sv.vcf |gzip -f > sv.gt.vcf.gz
# CNV calling
cnvpytor -root root.pytor -rd sort.cram -T GRCh38.fa && cnvpytor -root root.pytor -his 100 && cnvpytor -root root.pytor -partition 100 && cnvpytor -root root.pytor -call 100|perl cnv2vcf.pl -prefix id -fai GRCh38.fa.fai|bgzip -f > cnv.vcf.gz
# STR calling
ExpansionHunter --reads sort.cram --reference GRCh38.fa --variant-catalog variant_catalog.json --output-prefix ./kSTR 
```

# Linkag analysis

```
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar CombineGVCFs -R GRCh38.fa --variant A.gatk.gvcf.gz --variant B.gatk.gvcf.gz --variant C.gatk.gvcf.gz -O pedigree.gvcf.gz
java -Xmx3g -jar gatk-package-4.1.2.0-local.jar GenotypeGVCFs --QUIET true -R GRCh38.fa  --variant pedigree.gvcf.gz -O pedigree.vcf.gz
bcftools norm -m - pedigree.vcf.gz|vawk --header '{if($6>100 && I$AN==(NF-9)*2 && length($4)==1 && length($5)==1 && I$DP>(NF-9)*5 && I$MQ>55 && I$QD>5){print }}'|bgzip > pedigree.bi.vcf.gz

perl ./ensembl-vep-release-107/vep --assembly GRCh38 --fork 10 -i pedigree.bi.vcf.gz -o pedigree.vep.vcf.gz --vcf --compress_output bgzip --no_stats  --merged --force_overwrite --offline --use_given_ref \
--total_length --numbers --ccds --hgvs --symbol --canonical --protein --biotype --tsl --nearest symbol \
--fasta GRCh38.fa \
--dir_cache ./cache \
--dir_plugins ./VEP_plugins \
--custom gnomAD.vcf.gz,gnomADg,vcf,exact,0,AF,AF_eas,AF_popmax,popmax \
--custom mbiobank_ChinaMAP.phase1.vcf.gz,ChinaMAP,vcf,exact,0,AF,AC,AN

swift -p fam.txt --elod --penetrance=0.0,1.0,1.0 -f 0.0001 -u 10000 -o elod.txt 2>elod.txt
cat fam.txt|perl mogan.fam.pl > mogan.fam
tabix genetic.map.gz chr22|awk '{print $1,$2-1,$2}'|tr ' ' '\t'  > s.bed
tabix pedigree.vep.vcf.gz -R s.bed -h|perl filter.md.error.pl -f fam.txt|bcftools +split-vep - -s worst -f "%CHROM %POS %REF %ALT %ChinaMAP_AF %gnomADg_AF_eas[ %GT]\n"|perl mogan.ped.pl -minaf 0 -maxaf 0.01
plink --export ped --out tmp --tfile tmp
tabix genetic.map.gz -R mogan.s.bed|cut -f 4 > mogan.map
cat tmp.ped|cut -f '2,7-'|tr '\t' ' '|cat mogan.map tmp.af -|awk '{if(NR==1){print "map markers position";print ""};print}' > mogan.marker
gl_auto glautopar
tabix pedigree.vep.vcf.gz -R mogan.s.bed -h|bcftools +split-vep - -s worst -f "%CHROM %POS %REF %ALT %ChinaMAP_AF %gnomADg_AF_eas[ %GT]\n"|perl gigi.ped.pl -minaf 0 -maxaf 0.01
plink --export ped --out tmp --tfile tmp
cat tmp.ped|cut -f '2,7-'|tr '\t' ' ' > gigi.genotype
tabix genetic.map.gz -R gigi.s.bed|cut -f 4 > gigi.map
GIGI gigi.par
cat fam.txt|paste - impute.geno|tr '\t' ' '|cut -d ' ' -f '1-6,8-' > merlin.ped
tabix genetic.map.gz|awk '{if(NR==1){print "CHROMOSOME MARKER POSITION"};print $1" M"$2" "$4}' > merlin.map
tabix genetic.map.gz|awk '{if(NR==1){print "A Disease"};print "M M"$2}' > merlin.dat
Rscript lod.r merlin.ped merlin.map
```


