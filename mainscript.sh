##In twisst2 directory:
for chr in $(cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt); do
    TOPO=$(ls ${chr}| grep topocounts)
    INTERVAL=$(ls ${chr}| grep intervals)
    ##print max topo per row
    zcat ${chr}/${TOPO} | grep -v "#" | \
    awk 'BEGIN{ getline; for(i=1;i<=NF;i++) hdr[i]=$i; max=-1 }{ for(i=2; i<=NF; i+=2) if($i==max) { topo=topo","hdr[i] } else if($i>max) { max=$i; topo=hdr[i] } print topo; max=-1 }' \
    > ${chr}/maxtopo.tsv
    ## cursed one-liner: pastes intervals and topos together, uniq without sort to collapse consecutive ones; change end to realend (next start pos -1), remove "other" windows
    zcat ${chr}/${INTERVAL}| sed -n "2,$"p |paste - ${chr}/maxtopo.tsv |sed "s/,other//" | sed "s/other/topo80/" | uniq -f 3 |\
    awk 'BEGIN {FS=OFS="\t"}{chrom[NR] = $1; start[NR] = $2; end[NR] = $3; topo[NR] = $4} END {for (i = 1; i <= NR; i++) print chrom[i],start[i],end[i],topo[i],start[i+1]}'| \
    awk 'BEGIN {FS=OFS="\t"}{$6=$5-1;print $1,$2,$6,$4}'| sed -e "s/\.0//"| grep -v -E "," >${chr}/interval_topo.tsv &
done
wait

zcat chr29/Bos.*topocounts.tsv.gz |grep "#"| tr -d "#" >alltopos.tsv ##just finding newick trees of all the topologies from any random file
cat chr*/interval_topo.tsv > allintervals_complete_topo.tsv
cut -f 4 allintervals_complete_topo.tsv | sort |uniq -c| sort -nr > topofreq.txt
for topo in $(head -n 10 topofreq.txt|awk '{print $2}'); do t=$(grep "${topo} " alltopos.tsv|sed -E "s/\s/\t/"| tr -d ";"); f=$(grep "${topo}$" topofreq.txt | awk '{print $1}'|sed "s/\s//g"); echo -e "${t}\t${f}"; done > toptopos.tsv

awk 'BEGIN{FS=OFS="\t"}{if ($4=="topo34"||$4=="topo28"||$4=="topo88"||$4=="topo42"||$4=="topo20"||$4=="topo52"||$4=="topo36"||$4=="topo2"||$4=="topo54"){print $1,$2,$3,$4} else {print $1,$2,$3,"topo80"}}' allintervals_complete_topo.tsv |\
uniq -3 > allintervals_topo.tsv


###Create VCFs for each topology
mkdir topos
for topo in $(cut -f 1 toptopos.tsv); do
    mkdir topos/${topo}
    grep $topo allintervals_topo.tsv | cut -f 1,2,3 > topos/$topo/$topo.pos.txt
    mkdir topos/${topo}/vcfsubset
    for chr in $(cut -f 2 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt); do
        fullchr=$(grep "${chr}$" /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/reference/chromnames.txt| cut -f 1) ##real chrom name as indicated on ref genome (and in vcf)
        grep ${fullchr} topos/${topo}/${topo}.pos.txt > topos/${topo}/${chr}.${topo}.pos.txt ##find intervals for each chromosome for each topology by grepping real chrom name
        ##Extract snps at these intervals
        bcftools view --threads 8 -m2 -M2 -Oz -R topos/${topo}/${chr}.${topo}.pos.txt \
        /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/allmerged/Bos.${chr}.vcf.gz > topos/${topo}/vcfsubset/Bos.${chr}.vcf.gz &    
    done
done
wait

####snpcoords step--finds coordinates fulfilling various thresholds in each lineage
for topo in $(cut -f 1 toptopos.tsv); do
    cd topos/${topo}
    bash /home/export/base/ycsc_wx/liuz/online1/swl47/software/snpcoords.sh \
    "1.0" "0.3" "/home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/topos/${topo}/vcfsubset" \
    "/home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/lineage_nohybrids_modern.txt" \
    "/home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/slurmsub.txt"
    cd ../../
done

####snpfinder step--Traverses along tree from root to tip to determine lineages compared at each node, compares snpcoords output based on these to find diagnostic SNPs
#snpfinder
for topo in $(cut -f 1 toptopos.tsv); do
    nwk=$(grep "${topo}\s" toptopos.tsv | cut -f 2)
    cd topos/${topo}
    echo ${nwk} > Bos.${topo}.nwk
    sbatch /home/export/base/ycsc_wx/liuz/online1/swl47/software/snpfinder.sh "1.0" "0.3" "Bos.${topo}.nwk"
    cd ../../
done


##in aeDNA directory, where vcf from aeDNA is stored in ../vcf
####snpplacer step--assign lineages to samples based on pileup of SNPs at diagnostic sites
cd /home/export/base/ycsc_wx/liuz/online1/swl47/aeDNA/
mkdir twisst2; cd twisst2
for id in $(ls ../vcf | cut -f 1 -d . | sort -u ); do
    mkdir ${id}; cd ${id}
    for topo in $(cut -f 1 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/toptopos.tsv); do
        sbatch /home/export/base/ycsc_wx/liuz/online1/swl47/software/snpplacer.sh \
        "/home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/topos/${topo}" \
        "/home/export/base/ycsc_wx/liuz/online1/swl47/aeDNA/vcf/${id}.sorted.vcf.gz" ${topo}
    done
    cd ../
done

####Combine across output to get summary of assignments across all samples and all topologies
for id in $(ls ../vcf | cut -f 1 -d . | sort -u ); do
    echo "${id}" > ${id}/${id}_assignments.txt
    for topo in $(cut -f 1 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/toptopos.tsv); do
        node=$(tail -n 1 ${id}/${topo}/${topo}.pathfile.txt| cut -f 1| sed "s/_final//")
        snpsum=$(awk 'BEGIN{FS=OFS="\t"} {sum += $3} END {print sum}' ${id}/${topo}/${topo}.pathfile.txt)
        nextnode=$(grep "^${node}" /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/topos/${topo}/nodes.txt|cut -f 2)     
        ##find all lineages of nodes
        while [[ $(echo ${node}|grep [0-9]|wc -c) -gt 0 ]]; do
            for N in $(echo ${node}|grep -o "n[0-9]"); do
                match=$(grep "^$N" /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/topos/${topo}/nodes.txt| cut -f 2)
                nextnode2=$(echo ${nextnode}|sed "s/${N}/${match}/")
                nextnode=${nextnode2}
            done
            node=${nextnode}
        done
        echo "${node}_${snpsum}" >> ${id}/${id}_assignments.txt
    done
done
filenames=$(ls */*_assignments.txt | tr "\n" " " )
paste ${filenames} > allassignments.txt
echo "topologies" > alltopos.txt
for topo in $(cut -f 1 /home/export/base/ycsc_wx/liuz/online1/swl47/callsets/Bos/vcf/twisst2/toptopos.tsv);do echo ${topo} >>alltopos.txt; done
paste alltopos.txt allassignments.txt >alltoposallassignments.txt

