#!/bin/bash
mkdir scripts
sed -E -i "s/African [Tt]aurine/AfTau/" lineage_sticcs.txt
sed -i "s/India\-Pakistan indicine/SAsInd/" lineage_sticcs.txt #change to south asian indicine (they cluster tgt)
sed -i "s/North\-East Asian taurine/EAsTau/" lineage_sticcs.txt # change to East Asian indicine
sed -E -i "s/North European [Tt]aurine/EuTau/" lineage_sticcs.txt
sed -i "s/South-Central European taurine/EuTau/" lineage_sticcs.txt
sed -i -E "s/East Asian [Ii]ndicine/EAsInd/" lineage_sticcs.txt
sed -i "s/East Asian taurine/EAsTau/" lineage_sticcs.txt
sed -i "s/South Asian indicine/SAsInd/" lineage_sticcs.txt
cut -f 1 lineage_sticcs.txt > samples.txt
grep -v ERR4414056 lineage_sticcs.txt > lineage_twisst.txt
FILEPATH=../allmerged
for c in $(cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt | grep -E -v "X|Y|M")
do 
    mkdir ${c}
    cat slurmsub.txt > scripts/${c}_twisst2.sh
    echo "bcftools view --threads 60 -S samples.txt -o ${c}/Bos.${c}.sample.vcf.gz ${FILEPATH}/Bos.${c}.vcf.gz">>scripts/${c}_twisst2.sh
    echo "sticcs prep -i ${c}/Bos.${c}.sample.vcf.gz -o ${c}/Bos.${c}.sticcs.vcf.gz --outgroup ERR4414056">>scripts/${c}_twisst2.sh 
    echo "twisst2 -i ${c}/Bos.${c}.sticcs.vcf.gz -o ${c}/Bos.${c}.twisst2 --max_iterations 50 --ploidy 2 --group_names AfTau EAsInd EAsTau EuTau SAsInd --groups_file lineage_twisst.txt">>scripts/${c}_twisst2.sh  
    sbatch scripts/${c}_twisst2.sh 
done


c="chrM"
mkdir ${c}
cat slurmsub.txt > scripts/${c}_twisst2.sh
echo "bcftools view --threads 60 -S samples.txt -o ${c}/Bos.${c}.sample.vcf.gz ${FILEPATH}/Bos.${c}.vcf.gz">>scripts/${c}_twisst2.sh
echo "sticcs prep -i ${c}/Bos.${c}.sample.vcf.gz -o ${c}/Bos.${c}.sticcs.vcf.gz --outgroup ERR4414056">>scripts/${c}_twisst2.sh 
echo "twisst2 -i ${c}/Bos.${c}.sticcs.vcf.gz -o ${c}/Bos.${c}.twisst2 --max_iterations 50 --ploidy 1 --group_names AfTau EAsInd EAsTau EuTau SAsInd --groups_file lineage_twisst.txt">>scripts/${c}_twisst2.sh  
sbatch scripts/${c}_twisst2.sh 

##X chromosome--ploidy dependent on sex
#nano chrX_ploidy.txt
c="chrX"
mkdir ${c}
cat slurmsub.txt > scripts/${c}_twisst2.sh
echo "bcftools view --threads 60 -S samples.txt -o ${c}/Bos.${c}.sample.vcf.gz ${FILEPATH}/Bos.${c}.vcf.gz">>scripts/${c}_twisst2.sh
echo "sticcs prep -i ${c}/Bos.${c}.sample.vcf.gz -o ${c}/Bos.${c}.sticcs.vcf.gz --outgroup ERR4414056">>scripts/${c}_twisst2.sh 
echo "twisst2 -i ${c}/Bos.${c}.sticcs.vcf.gz -o ${c}/Bos.${c}.twisst2 --max_iterations 50 --ploidy_file chrX_ploidy.txt --group_names AfTau EAsInd EAsTau EuTau SAsInd --groups_file lineage_twisst.txt">>scripts/${c}_twisst2.sh  
sbatch scripts/${c}_twisst2.sh 

##Y chromosome--only males, haploid
##nano lineage_twisst.males.txt
##cut -f 1 lineage_twisst.males.txt > samples.males.txt
c="chrY"
mkdir ${c}
cat slurmsub.txt > scripts/${c}_twisst2.sh
echo "bcftools view --threads 60 -S samples.males.txt -o ${c}/Bos.${c}.sample.vcf.gz ${FILEPATH}/Bos.${c}.vcf.gz">>scripts/${c}_twisst2.sh
echo "sticcs prep -i ${c}/Bos.${c}.sample.vcf.gz -o ${c}/Bos.${c}.sticcs.vcf.gz --outgroup ERR4414056">>scripts/${c}_twisst2.sh 
echo "twisst2 -i ${c}/Bos.${c}.sticcs.vcf.gz -o ${c}/Bos.${c}.twisst2 --max_iterations 50 --ploidy 1 --group_names AfTau EAsInd EAsTau EuTau SAsInd --groups_file lineage_twisst.males.txt">>scripts/${c}_twisst2.sh  
sbatch scripts/${c}_twisst2.sh 
