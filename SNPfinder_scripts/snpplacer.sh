#!/bin/bash
#SBATCH --nodes=2
#SBATCH --cpus-per-task=60
#SBATCH --time=2:00:00
#SBATCH --mem=20GB
#SBATCH --mail-user=swl47@cam.ac.uk
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --job-name=snpplacer
#SBATCH --error=tmp_snpplacer/snpplacer.%J.err
#SBATCH --output=tmp_snpplacer/snpplacer.%J.out

SNPFINDER=$1
VCF=$2
ID=$3

node=$(tail -n 1 ${SNPFINDER}/nodes.txt | cut -f 1)
mkdir ${ID}
confi=1
echo -e "${node}\tNA" >${ID}/${ID}.pathfile.txt
while [[ $(echo $node | wc -c) == 3 ]]; do ##no of characters =3, n1 has 3 characters as wc -c counts newline
    g1=$(grep -E "^${node}" ${SNPFINDER}/nodes.txt|cut -f 2|cut -f 1 -d ,) ##group 1
    g2=$(grep -E "^${node}" ${SNPFINDER}/nodes.txt|cut -f 2|cut -f 2 -d ,) ##group 2
    cut -f 1,2 ${SNPFINDER}/snps/${node}.snps.txt |grep -E 'NC_[0-9]+\.1'$'\t[0-9]{1,9}$'> ${ID}/${node}.pos.txt #extract positions for node , make sure formatting is right
    bcftools query -R ${ID}/${node}.pos.txt -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]' ${VCF} >> ${ID}/${ID}.${node}.txt 
    cut -f 1,2,3 ${ID}/${ID}.${node}.txt | xargs -I "?" grep "?" ${SNPFINDER}/snps/${node}.snps.txt > ${ID}/${ID}.${node}.nodesnps.txt
    ##paste snp file and sample GTs together, (they should be in same order), filter to make sure ref and alt are the same, print chr pos ref alt n1 n2 genotype; then assign allele based on genotype, then check for g1 or g2
    paste ${ID}/${ID}.${node}.nodesnps.txt ${ID}/${ID}.${node}.txt | awk 'BEGIN {FS=OFS="\t"}{if ($1==$7 && $2==$8 && $3==$9 && $4==$10) {print $1,$2,$3,$4,$5,$6,$11}}' |\
    awk 'BEGIN {FS=OFS="\t"}{if ($7=="0/0") {print $1,$2,$3,$4,$5,$6,$3} else if ($7=="1/1") {print $1,$2,$3,$4,$5,$6,$4}}'|\
    awk 'BEGIN {FS=OFS="\t"}{if ($7==$5) {print $1,$2,$3,$4,$5,$6,$7,"g1"} else if ($7==$6) {print $1,$2,$3,$4,$5,$6,$7,"g2"} else {print $1,$2,$3,$4,$5,$6,$7,"node"}}'\
    >${ID}/${ID}.${node}.filtered.txt
    ##number of group 1 and number of group2
    nog1=$(cut -f 8 ${ID}/${ID}.${node}.filtered.txt | grep g1 | wc -l )
    nog2=$(cut -f 8 ${ID}/${ID}.${node}.filtered.txt | grep g2 | wc -l )
    ##compare number of group 1 vs number of group 2, print the one that has more
    if [[ ${nog1} -gt ${nog2} ]]; then node=$g1; weight=$(echo "scale=5 ; $nog1 / ($nog1 + $nog2)" | bc ); winno=${nog1};
    elif [[ ${nog2} -gt ${nog1} ]]; then node=$g2; weight=$(echo "scale=5 ; $nog2 / ($nog1 + $nog2)" | bc); winno=${nog2};
    else nfinal=$(echo ${node}_final); node=$nfinal; weight="0.5"
    fi
    echo -e "${node}\t${weight}\t${winno}" >> ${ID}/${ID}.pathfile.txt
done
