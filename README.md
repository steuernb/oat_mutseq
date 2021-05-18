# Mutational Genomics to identify candidate genes in oat

We used these scripts to first re-clone SAD3 as proof of concept, then use the method to find a candidate for SAD4-like (PAL2).

The reference genome sequence can be found here: [http://db.ncgr.ac.cn/oat/db/data2/OAT_v2_final_chr.fa.gz](http://db.ncgr.ac.cn/oat/db/data2/OAT_v2_final_chr.fa.gz)

The rawdata of re-sequenced mutants is available from European Nucleotide Archiv, study number [PRJEB45039](http://www.ebi.ac.uk/ena/data/view/PRJEB45039). (publication date: 18-Jul-2021)




## Method



### Read mapping

#### Prerequisites

BWA: [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)
SAMTOOLS: [http://www.htslib.org/](http://www.htslib.org/)

#### Usage

Index the reference

```
bwa index OAT_v2_final_chr.fa
samtools faidx OAT_v2_final_chr.fa
```

Run the mappings for each read pair.

```
bwa mem OAT_v2_final_chr.fa mutant1_read1.fq.gz mutant1_read2.fq.gz | samtools sort -o mutant1.bam -
```

Some samples have more than one read pair. Run mappings individually, then merge resulting bam files

```
samtools merge mutant1.bam mutant1_sub1.bam mutant1_sub2.bam ...
```

Index BAM files and run mpileup

```
samtools index -c mutant1.bam
samtools mpileup -BQ0 -f OAT_v2_final_chr.fa mutant1.bam > mutant1.mpileup
```

### Call SNPs

Run [QuickSNP](https://github.com/steuernb/oat_mutseq/blob/main/QuickSNP.jar) from this repository to pre-select mutated positions. Do this for each mutant.

```
java -jar QuickSNP.jar -i mutant1.mpileup -o mutant1.snp.txt
```

### Select candidates

Use HuntSNPs with various parameter settings to identify candidate genes. For each potential candidate, manually check mappings using a genome browser.

We found the easiest way is to start varying the coverage parameter, starting from a very high threshold to a lower. In our case, this was sufficient to identify the candidate, however, depending on how noisy data is, allele frequency stringency could be reduced or the number of mutants necessary to report a region can be lowered.

```
java -jar HuntSNPs.jar -i mutant1.snp.txt mutant2.snp.txt mutant3.snp.txt mutant4.snp.txt -w 10000 -n 4 -e 4 -a 0.1 -s 1 -c 20
```

Parameter | Description
--- | ---
-i | Input files. These are output files of QuickSNP
-w | Window size. Size of interval to find mutations (default 10000)
-n | Minimum number of mutants with mutation in an interval (default: 4)
-e | Minimum number of mutants with canonical Sodium Azide mutation in an interval (default: 4)
-a | Maxium reference allele frequency to consider a SNP (default: 0.1)
-c | Minimum coverage to consider a SNP (default: 20)
-s | Maximum mutants allowed to share a position (default: 1)
					






