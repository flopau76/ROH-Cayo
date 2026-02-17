# Code for Inferring ROH on Cayo Santiago
### Florence Pautet

## Preparing the vcf with the chinese samples
### Downloading the data
```bash
# Get the vcf from gigaDB
wget https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100484/population.vcf.gz -O 81wildChineseRhesus.vcf.gz
```
### Filtering the data
```bash
# 1) Remove CR1 and CR2 individuals & Filter multi-allelic snps and indels; keep only passed sites
# 2) Set GT to missing with GQ<20 or lowDP flag
# 3) Filter missing GT
# 4) Rename chromosomes from 'X' to 'chrX' for compatibility with the chain file
# 5) Recompress the output
bcftools view -s ^CR1,CR2 81wildChineseRhesus.vcf.gz -m2 -M2 -v snps -i 'FILTER="PASS"' -Ou |
bcftools filter -i 'FMT/GQ>19 & FMT/FT!="lowDP"' -S . -Ou |
bcftools view -i 'F_MISSING<0.2' -Ov |
awk 'BEGIN { OFS="\t" }
     /^##contig=<ID=[^,]+/ {
          sub(/^##contig=<ID=/, "##contig=<ID=chr"); print; next
     }
     /^#/ { print; next }
     {$1="chr"$1; print}' |
bcftools view -Oz -o 79wildChineseRhesus.vcf.gz
```
### Lifting over to Mmul10
#### Preparing the files
```bash
# Download the reference genome and the chain file
wget https://hgdownload.gi.ucsc.edu/goldenPath/rheMac8/liftOver/rheMac8ToRheMac10.over.chain.gz
wget https://hgdownload.gi.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz
# Change all softmasking in the reference to uppercase for compatibility with the chain file
zcat rheMac10.fa.gz | awk '/^>/ {print; next} {print toupper($0)}' | bgzip > rheMac10_upper.fa.gz
# Create sequence dictionary for the reference genome
gatk CreateSequenceDictionary -R rheMac10_upper.fa.gz -O rheMac10_upper.dict
```

#### Lifting over
```bash
# Lift over
gatk LiftoverVcf \
--INPUT 79wildChineseRhesus.vcf.gz \
--OUTPUT 79wildChineseRhesus_mmul10.vcf.gz \
--CHAIN rheMac8ToRheMac10.over.chain.gz \
--REJECT rejected_variants.vcf.gz \
--REFERENCE_SEQUENCE rheMac10_upper.fa.gz

# Filter out unmapped contigs and rename chromosomes from 'chrX' to 'X'
(for i in {1..20}; do echo "chr$i $i"; done;) > mapping.txt
bcftools view -r $(for i in {1..20}; do echo -n "chr$i,"; done;) -Ou 79wildChineseRhesus_mmul10.vcf.gz |
bcftools annotate --rename-chrs mapping.txt -Oz -o 79wildChineseRhesus_mmul10_final.vcf.gz
rm mapping.txt
```

## Preparing the vcf with the Cayo Santiago samples
TODO

## Calling the roh
### Applying bcftools/roh
```bash
# Call on chinese genomes
bcftools roh 79wildChineseRhesus_mmul10_final.vcf.gz --AF-tag AF -m genetic_map/chr{CHROM}.txt --output-type r --output 79wildChineseRhesus.rg
# Call on Cayo Santiago genomes
bcftools roh 98CayoRhesus.vcf.gz --AF-tag AF -m genetic_map/chr{CHROM}.txt --output-type r --output 98CayoRhesus.rg
```
### Postprocessing the files
```bash
python3 post_process.py -i 79wildChineseRhesus.rg -m genetic_map/chr{CHROM}.txt
python3 post_process.py -i 98CayoRhesus.rg -m genetic_map/chr{CHROM}.txt
```

```bash
```