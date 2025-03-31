#!/bin/bash
#SBATCH --nodes=2
#SBATCH --cpus-per-task=60
#SBATCH --time=1-02:00:00
#SBATCH --mem=20GB
#SBATCH --mail-user=swl47@cam.ac.uk
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --job-name=snpfinder
#SBATCH --error=tmp_snpfinder/snpfinder.%J.err
#SBATCH --output=tmp_snpfinder/snpfinder.%J.out
THRESHOLD=$1
t=$2
TREEFILE=$3
#th=$(echo $THRESHOLD| cut -f 2 -d .)
#t=$(echo $((10-th))| sed -e "s/^/0\./")
##Find nodes
cat $TREEFILE > trees.txt
i=1
N=$(grep -o , $TREEFILE | wc -l) ##number of nodes
lasttree=$(tail -n 1 $TREEFILE) ##latest tree
echo -e "nodes\tlineages" > nodes.txt
while [[ $i -le $N ]]
 do
    for node in $(tail -n1 trees.txt | grep -E -o "([A-Za-z0-9]+,[A-Za-z0-9]+)" ); do 
        echo -e "n${i}\t$node" >> nodes.txt
        lasttree_tmp=$(echo $lasttree|sed "s/(${node})/n${i}/"); lasttree=$(echo $lasttree_tmp); echo $lasttree>> trees.txt; ((i+=1)); done
done

##find diagnostic SNPs at each node
mkdir snps
##make node coords; then find common
for line in $(seq 2 $(wc -l nodes.txt|cut -f 1 -d " ")); do ## from seq 2 due to header
    node=$(sed -n "${line}p" nodes.txt | cut -f 1)
    g1=$(sed -n "${line}p" nodes.txt | cut -f 2| cut -f 1 -d ,)
    g2=$(sed -n "${line}p" nodes.txt | cut -f 2| cut -f 2 -d ,)
    #echo -e "chr\tpos\tREF\tALT\t$g1\t$g2" > snps/$node.snps.txt
    for chrcoord in $(ls snpcoords/| grep coords); do
        comm -1 -2 snpcoords/$chrcoord/$g1.${THRESHOLD}.alt.coord snpcoords/$chrcoord/$g2.${THRESHOLD}.alt.coord > snpcoords/$chrcoord/$node.${THRESHOLD}.alt.coord &
        pids1[${chrcoord}]=$! 
        comm -1 -2 snpcoords/$chrcoord/$g1.${t}.alt.coord snpcoords/$chrcoord/$g2.${t}.alt.coord > snpcoords/$chrcoord/$node.${t}.alt.coord &
        pids2[${chrcoord}]=$! 
        comm -1 -2 snpcoords/$chrcoord/$g1.${THRESHOLD}.ref.coord snpcoords/$chrcoord/$g2.${THRESHOLD}.ref.coord > snpcoords/$chrcoord/$node.${THRESHOLD}.ref.coord &
        pids3[${chrcoord}]=$!    
        comm -1 -2 snpcoords/$chrcoord/$g1.${t}.ref.coord snpcoords/$chrcoord/$g2.${t}.ref.coord > snpcoords/$chrcoord/$node.${t}.ref.coord &
        pids4[${chrcoord}]=$!
    done
    for pid1 in ${pids1[*]}; do wait $pid1; done
    for pid2 in ${pids2[*]}; do wait $pid2; done
    for pid3 in ${pids3[*]}; do wait $pid3; done
    for pid4 in ${pids4[*]}; do wait $pid4; done
    #g1 alt, g2 ref
    ls snpcoords|grep coords| xargs -I ? comm -1 -2 snpcoords/"?"/$g1.${t}.alt.coord snpcoords/"?"/$g2.${THRESHOLD}.ref.coord | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$4,$3}' >> snps/$node.snps.txt & ##g1 is alt, g2 is ref format is chr pos ref alt g1(alt) g2(ref)
    ls snpcoords|grep coords| xargs -I ? comm -1 -2 snpcoords/"?"/$g1.${THRESHOLD}.alt.coord snpcoords/"?"/$g2.${t}.ref.coord | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$4,$3}' >> snps/$node.snps.txt & ##g1 is alt, g2 is ref format is chr pos ref alt g1(alt) g2(ref)
    #g1 ref, g2 alt
    ls snpcoords|grep coords| xargs -I ? comm -1 -2 snpcoords/"?"/$g1.${t}.ref.coord snpcoords/"?"/$g2.${THRESHOLD}.alt.coord | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$3,$4}' >> snps/$node.snps.txt & ##g1 is ref, g2 is alt format is chr pos ref alt g1(ref) g2(al)
    ls snpcoords|grep coords| xargs -I ? comm -1 -2 snpcoords/"?"/$g1.${THRESHOLD}.ref.coord snpcoords/"?"/$g2.${t}.alt.coord | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$3,$4}' >> snps/$node.snps.txt & ##g1 is ref, g2 is alt format is chr pos ref alt g1(ref) g2(al)
done
wait
