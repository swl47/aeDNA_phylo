#!/bin/bash
##all
cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt|\
xargs -I ? plink2 --vcf allmerged/Bos."?".vcf.gz --keep allid.txt --allow-extra-chr --mind --min-alleles 2 --max-alleles 2 --rm-dup exclude-all --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out pca/Bos."?"
##prune and create pca
cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt|\
xargs -I ? plink2 --vcf allmerged/Bos."?".vcf.gz --keep allid.txt  --allow-extra-chr --mind --min-alleles 2 --max-alleles 2 --set-missing-var-ids @:# --extract pca/Bos."?".prune.in --make-bed --pca --out pca/Bos."?"

##modern, non-hybrids
cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt|\
xargs -I ? plink2 --vcf ../allmerged/Bos."?".vcf.gz --keep nohybrids_modern_id.txt --allow-extra-chr --mind --min-alleles 2 --max-alleles 2 --rm-dup exclude-all --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Bos."?"
##prune and create pca
cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt|\
xargs -I ? plink2 --vcf ../allmerged/Bos."?".vcf.gz --keep nohybrids_modern_id.txt --allow-extra-chr --mind --min-alleles 2 --max-alleles 2 --set-missing-var-ids @:# --extract Bos."?".prune.in --make-bed --pca --out Bos."?"
