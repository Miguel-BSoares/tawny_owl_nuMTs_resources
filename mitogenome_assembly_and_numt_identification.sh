####Pipeline to extract, parse and assemble mitochondrial genome##### ALWAYS SAME INDIVIDUAL , NO CROSS MAPPING
#blast pacbio reads against multiple mitogenomes
for file in to_blast/*.fasta; do blasr /scratch/project_2004038/owl_data/m64048_210523_061850.subreads.bam $file -m 1 > $file.out;
done
blasr owl_data/sample.bam seed0.fasta -m 1 --hitPolicy all > blast_data/v6seed1.blasr.out
#create fastq file with reads to assemble - mito_seqkt.sh
seqtk subseq /scratch/project_2004038/owl_data/pacbio.V5A.subreads.fastq.gz v5seed2.blasr.sort.list.out.fixed > v5mito.fastq
#multiple canu assemblies - 10 here - canu_loop.sh
mkdir canu_assemblies_v5
for (( i = 1; i<= 10;i++ ))
do
	mkdir -p mito_canu$i
	canu -d mito_canu$i -p owl genomeSize=16k -pacbio v5mito.fastq correctedErrorRate=0.035 utgOvlErrorRate=0.065 trimReadsCoverage=2 trimReadsOverlap=500 MhapSensitivity=low purgeOverlaps=aggressive useGrid=remote gridEngineArrayMaxJobs=20
	mv mito_canu$i/owl.contigs.fasta canu_assemblies_v5/assembly$i.fasta
done

######Trycyler - software to make consensus out of different assemblies and circularize genomes; different steps ##### 
trycycler cluster --assemblies canu_assemblies_v5/*.fasta --reads v5mito.fastq --out_dir trycycler_output_v5
##choose on which clusters to progress with to the next step###
##help is provided by the newick tree##
##eliminate duplicated 
##trycycler_reconcile.sh
trycycler reconcile --reads v6mito.fastq --cluster_dir /scratch/project_2004038/blast_data/1_contigs/cluster_001 --verbose
##trycycler_msa.sh
trycycler msa --cluster_dir /scratch/project_2004038/blast_data/1_contigs/cluster_001
##trycycler_parti.sh
trycycler partition --reads v6mito.fastq --cluster_dir /scratch/project_2004038/blast_data/1_contigs/cluster_001
##trycycler_cons.sh
trycycler consensus --cluster_dir /scratch/project_2004038/blast_data/1_contigs/cluster_001
## polishing o trycycler final assembly - pmmb2_mito.sh & gcpp_mito.sh
export PROJAPPL=/projappl/project_2004038
module load bioconda
conda activate pacbio
mkdir -p pbmm2_mito
export TMPDIR=/scratch/project_2004038/pbmm2_mito/
pbmm2 index blast_data/1_contigs/cluster_001/7_final_consensus.fasta pbmm2_mito/mito_assembly.mmi
pbmm2 align --sort --log-level DEBUG pbmm2_mito/mito_assembly.mmi owl_data/m64048_210523_061850.subreads.bam pbmm2_mito/mito.bam
mkdir -p mito_polished
gcpp pbmm2_mito/mito.bam -r blast_data/1_contigs/cluster_001/7_final_consensus.fasta -o mito_polished/mito_polished.fasta --algorithm=arrow -j 10
#align mitogenomes and organize results file to compare score and % identity - all2_blasr.sh
module load gcc/9.1.0
module load blast/2.12.0
mtgenomes="Saluc_ncbi
saluc_pacbio_v5
saluc_pacbio_v6
Slepto
Socci
Suralen
Svaria
Talba
Tlong"
for sample in $mtgenomes; do blasr *.fasta ${sample}.fasta -m 1 | sort | awk '{print $1,$2,$5,$6}' >>all.blasr.out
done

#####capture reads, from raw pacbio data, 
#blast pacbio reads against multiple mitogenomes and keeping frag lenghts between 100bp and 20000bp (approx. max mitogenome size) - 
module load gcc/9.1.0
module load blast/2.12.0
for file in to_blast/*.fasta; do blasr /scratch/project_2004038/owl_data/m64048_210523_061850.subreads.bam $file -m 1 --hitPolicy all > $file.out;
done
##sort hits by size 
for file in  ./*out;do awk '{a=$8-$7;print $0,a;}' $file | sort -n -r -k14,14 | awk '$14>100 && $14<20000' | sort -uk1,1 > $file.sort.out;
done
#retrieve data on alignement score and alignment % identity
for file in ./*.sort.out; do awk '{print $1,$5,$6}' $file > $file.align.stats; 
done

###CROSS MAPPING STRATEGY BETWEEN RAW PACBIO READS AND CLEANED MITOGENOMES#####can be incorporated on the previous step after both circular genomes are available for grey and brown morphs
##utilize blasr to fish for subreads that match the mitogenome
blasr subreads_file mitogenome_file.gz -m 1 --hitPolicy all > output
#utilize seqtk to extract sequences based on a file with coordinate information
seqtk subseq {GENOME_FILE}.gz coordinates_file > output.fastq


