Sorry, this is late but I realised I should add a README so here it is:  
  
Files in the Snakemake folder are Snakefiles that run the variant calling pipeline  
Snakemake_panel, followed by the merge scripts, was used for the reference panel and calls from raw FASTQ reads. Snakemake_aeDNA was the pipeline for calling variants from the aeDNA data (which was provided to me as already-aligned bam files)  
  
Files in the SNPfinder_scripts folder are scripts used in various steps of the SNPfinder pipeline  
mainscript.sh calls the scripts in SNPfinder_scripts to run the whole SNPfinder pipeline, from ARG inference and topology weighting, to final lineage assignments  
  
...My apologies for writing most of this in bash, I regret that decision too
