#!/bin/bash
THRESHOLD=$1
t=$2
FILEPATH=$3
LINEAGE=$4
SLURMSUB=$5
mkdir snpcoords
mkdir snpcoords/tmp
#th=$(echo $THRESHOLD| cut -f 2 -d .)
 ##THRESHOLD=0.8, t =0.2
#t=$(echo $((10-th))| sed -e "s/^/0\./")
for chr in $(cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt); do
    file="${FILEPATH}/Bos.${chr}.vcf.gz"
    mkdir snpcoords/${chr}_coords
    for g in $(cut -f 2 $LINEAGE| sort | uniq)
    do
        grep ${g} $LINEAGE | cut -f 1 > ${g}_samples.txt ###SRRs for given group
        cat ${SLURMSUB} > snpcoords/tmp/${chr}.${g}.coord.sh
        echo "bcftools view --threads 8 -m2 -M2 -S ${g}_samples.txt -Ou $file | bcftools query -e 'GT =\"AA\" || GT =\"AR\" || GT =\"RA\"' -f '%CHROM\t%POS\t%REF\t%ALT\n' - | sort > snpcoords/${chr}_coords/${g}.${THRESHOLD}.ref.coord" >> snpcoords/tmp/${chr}.${g}.coord.sh 
        echo "bcftools view --threads 8 -m2 -M2 -S ${g}_samples.txt -Ou $file | bcftools query -i 'F_PASS(GT=\"RR\") > $t' -f '%CHROM\t%POS\t%REF\t%ALT\n' - | sort > snpcoords/${chr}_coords/${g}.${t}.ref.coord" >> snpcoords/tmp/${chr}.${g}.coord.sh 
        echo "bcftools view --threads 8 -m2 -M2 -S ${g}_samples.txt -Ou $file | bcftools query -e 'GT =\"RR\" || GT =\"AR\" || GT =\"RA\"' -f '%CHROM\t%POS\t%REF\t%ALT\n' - | sort > snpcoords/${chr}_coords/${g}.${THRESHOLD}.alt.coord" >> snpcoords/tmp/${chr}.${g}.coord.sh 
        echo "bcftools view --threads 8 -m2 -M2 -S ${g}_samples.txt -Ou $file | bcftools query -i 'F_PASS(GT=\"AA\") > $t' -f '%CHROM\t%POS\t%REF\t%ALT\n' - | sort > snpcoords/${chr}_coords/${g}.${t}.alt.coord" >> snpcoords/tmp/${chr}.${g}.coord.sh 
        sbatch snpcoords/tmp/${chr}.${g}.coord.sh
    done
done
