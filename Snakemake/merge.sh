##Scripts for merging vcf files; first one merges across one accession, second one combines the files acrosss all accessions

#!/bin/bash
##mergevcf.sh; merge vcfs for each accession
#bash mergevcf.sh "name"
mkdir mergevcf
mkdir mergevcf/tmp
ls vcf/*.sorted.vcf.gz| grep gz > vcfacc.txt
file=$(ls vcf/*.sorted.vcf.gz|grep gz|sed -n 1p)
for chr in $(bcftools index -s ${file}| cut -f 1 | grep "NC"); do chrname=$(grep ${chr} /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt | cut -f 2);
cat /home/export/base/ycsc_wx/liuz/online1/swl47/software/mergevcf_slurmsub.txt > mergevcf/tmp/merge${chrname}.sh
echo "bcftools merge --threads 60 -l vcfacc.txt -r ${chr} | bcftools filter -W --threads 60 -i'QUAL>20 && DP>10' -o mergevcf/${1}.${chrname}.filt.vcf.gz">> mergevcf/tmp/merge${chrname}.sh
sbatch mergevcf/tmp/merge${chrname}.sh
done

##mergeall.sh
#!/bin/bash
mkdir tmp
chrfile=/home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt
for chr in $(cat ${chrfile}| cut -f 1); do 
    chrname=$(grep ${chr} ${chrfile} | cut -f 2)
    ls sep_vcfs/*.${chrname}.filt.vcf.gz > tmp/${chrname}.vcf.txt
    cat slurmsub.txt > tmp/merge${chrname}.sh
    echo "bcftools merge --threads 60 -l tmp/${chrname}.vcf.txt | bcftools view -i 'F_MISSING<0.3 && MAF>0.05' -W --threads 60 -o allmerged_19_2/Bos.${chrname}.vcf.gz" >> tmp/merge${chrname}.sh
    sbatch tmp/merge${chrname}.sh
done
