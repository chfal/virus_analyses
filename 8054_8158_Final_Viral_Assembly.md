# 8054 and 8158 Final Assembly

Here is the code associated for the paper "Whole-Genome Sequencing and Phylogenetic Analysis of Anolis Adenovirus 2 Reveals Conserved Genome Organization and Gene-Specific Evolutionary Patterns"; Falvey and Geneva; 2026. In total, we did 3 rounds of sequencing for sample 8054 and 2 rounds of sequencing for sample 8158, which were done to make the coverages more even.

# Preparing Reads for Assembly

## Sample 8054 and 8158 FastQC

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=fastqc                       # job name for listing in queue
#SBATCH --mem=50G                              # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module load java                                # load any modules needed
module load FastQC

echo "Bash commands for the analysis you are going to run"

readset="/projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001"
readset="/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001"
readset="/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001"
echo "##################### fastqc initial quality analysis"
fastqc -t 20 \
${readset}_R1_001.fastq \
${readset}_R2_001.fastq \
-o /projects/f_geneva_1/chfal/virus_analyses/fastqc/run_two
```

## Trimmomatic to trim low-quality reads for 8054 and 8158

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=fastqc                       # job name for listing in queue
#SBATCH --mem=20G                              # memory to allocate in Mb
#SBATCH -n 1                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module load java                                # load any modules needed
module load FastQC

echo "Bash commands for the analysis you are going to run"

#readset="/projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001"
#readset="/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001"
readset="/projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001"

echo "##################### fastqc initial quality analysis"
echo "##################### trimmomatic"
java -jar /projects/f_geneva_1/programs/trimmomatic/trimmomatic-0.39.jar PE \
-threads 20 -phred33 -trimlog ${readset}_trim.log \
${readset}_R1_001.fastq ${readset}_R2_001.fastq \
${readset}_filtered.R1.fq ${readset}_filtered.unpaired.R1.fq \
${readset}_filtered.R2.fq ${readset}_filtered.unpaired.R2.fq \
ILLUMINACLIP:/projects/f_geneva_1/programs/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:4 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:13:20 MINLEN:23



echo ""
echo "##################### fastqc trimmomatic quality analysis"
fastqc -t 20 \
${readset}_filtered.R1.fq \
${readset}_filtered.R2.fq \
-o /projects/f_geneva_1/chfal/virus_analyses/8054_2026

```

# Genome Assembly

We wanted to do four assemblies total, 8054/8158 in both careful and metaviral mode. It was also important that everything got output to scratch as Spades in careful mode will actually completely crash the cluster if it is going to a projects folder as it will make more than 1.5 million files (it makes an error correction file for every single contig; there are sometimes more than 700,00 contigs).

It's important to note the headers of these slurm jobs, as some were not making it all the way through the assembly process with 250GB and 3 days. Sometimes it just goes and decides to crash which is odd/annoying. It should take less than 24 hours; sometimes it doesn't finish. If it does crash or stop working, there is no point in checkpointing it, the documentation will say you can, but it doesn't actually work, so just move on and resubmit this.


## 8054 Careful

```
#!/bin/bash
#SBATCH --partition=p_geneva_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=spades                       # job name for listing in queue
#SBATCH --mem=250G                              # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=14-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


module purge
eval "$(conda shell.bash hook)"
conda activate SPAdes
spades.py --careful -t 20 -1 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R1.fq -1 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R1.fq -1 /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R1.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R2.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R2.fq -2 /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R2.fq  -o /scratch/chf29/8054_final
```

## 8158 Careful

```
#!/bin/bash
#SBATCH --partition=mem                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002,halc068               # exclude CCIB GPUs
#SBATCH --job-name=spades                       # job name for listing in queue
#SBATCH --mem=250G                              # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


module purge
eval "$(conda shell.bash hook)"
conda activate SPAdes
# spades.py --careful -t 20 -1 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R1.fq -1 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R1.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R2.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R2.fq -o /scratch/chf29/8158_spades

spades.py --meta -t 20 -1 /scratch/chf29/8158_S2_L001_filtered.R1.fq -1 /scratch/chf29/8158_S2_L001_filtered.R1_v2.fq -2 /scratch/chf29/8158_S2_L001_filtered.R2.fq -2 /scratch/chf29/8158_S2_L001_filtered.R2_v2.fq -o /scratch/chf29/8158_meta
```



## 8054 Metaviral 

```
#!/bin/bash
#SBATCH --partition=mem                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=spades                       # job name for listing in queue
#SBATCH --mem=250G                              # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


module purge
eval "$(conda shell.bash hook)"
conda activate SPAdes
spades.py --meta -t 20 -1 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R1.fq -1 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R1.fq -1 /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R1.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R2.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R2.fq -2 /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R2.fq  -o /scratch/chf29/8054_meta
```

## 8158 Metaviral 

```

#!/bin/bash
#SBATCH --partition=mem                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002,halc068               # exclude CCIB GPUs
#SBATCH --job-name=spades                       # job name for listing in queue
#SBATCH --mem=250G                              # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE      # email for the following reasons


module purge
eval "$(conda shell.bash hook)"
conda activate SPAdes
# spades.py --careful -t 20 -1 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R1.fq -1 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R1.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R2.fq -2 /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R2.fq -o /scratch/chf29/8158_spades

spades.py --meta -t 20 -1 /scratch/chf29/8158_S2_L001_filtered.R1.fq -1 /scratch/chf29/8158_S2_L001_filtered.R1_v2.fq -2 /scratch/chf29/8158_S2_L001_filtered.R2.fq -2 /scratch/chf29/8158_S2_L001_filtered.R2_v2.fq -o /scratch/chf29/8158_meta

```

# Blast

We then had to blast all the contigs out. This is a fairly standard Blast command so I just changed the file paths as necessary.

```
#!/bin/bash
#SBATCH --partition=main                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=blast                       # job name for listing in queue
#SBATCH --mem=30G                              # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=FAIL      # email for the following reasons


module purge
module load blast/2.10.1-zz109
blastn -query /projects/f_geneva_1/chfal/virus_analyses/8054_final/scaffolds.fasta -db /projects/ccib/shain/blastdb/nt -evalue 1e-10 -num_threads 32 -out /projects/f_geneva_1/chfal/virus_analyses/8054_final/8054_final_blast_results.txt -max_target_seqs 5 -outfmt "6 sscinames qseqid sseqid pident length mismatch evalue bitscore"
```

In all of the Spades output, the largest final scaffold was around ~37KB. This is perfect size of the virus. And no surprise, in each assembly the largest contig matched as Anolis adenovirus / Lizard Adenovirus. It is clear that we assembled the largest contig and the virus almost perfectly. There are also some smaller contigs (~4000bp) that also are BLASTING as positive to the virus but if I recall we just select the largest contig and use that as our genome assembly as these viral fragments are not able to be placed successfully.

Example code to get largest scaffold and its contents:

```
awk '/^>/{if(N++)exit} {print}' scaffolds.fasta > 8054_careful_first_scaffold.fasta
```

Once we BLAST and find the largest scaffold thankfully happens to be Anolis Atadenovirus, it's time to do BWA to check the assembly coverage and depth to make sure we don't need to do more sequencing. As such:


# BWA and Coverage Depth

## 8054 careful Bam

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=bwa                          # job name for listing in queue
#SBATCH --mem=50G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

echo "load any Amarel modules that script requires"
module purge # clears out any pre-existing modules
module load samtools # load any modules needed
module load bwa


ROUND_1=/projects/f_geneva_1/chfal/virus_analyses/aug_23_sequencing/round_1
ROUND_2=/projects/f_geneva_1/chfal/virus_analyses/aug_23_sequencing/round_2
# FILE=$1

echo "Bash commands for the analysis you are going to run"

echo "##################### index and align with BWA"
bwa index /projects/f_geneva_1/chfal/virus_analyses/8054_final/8054_careful_first_scaffold.fasta

bwa mem -t 10 -5SPM /projects/f_geneva_1/chfal/virus_analyses/8054_final/8054_careful_first_scaffold.fasta \
<(cat /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R1.fq /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R1.fq /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R1.fq) \
<(cat /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R2.fq /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R2.fq /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R2.fq) | \
samtools sort -@10 -o /projects/f_geneva_1/chfal/virus_analyses/8054_final/mapped_8054_careful.bam -
samtools faidx /projects/f_geneva_1/chfal/virus_analyses/8054_final/8054_careful_first_scaffold.fasta

echo "change user group of files created"
chgrp -R g_geneva_1 /projects/f_geneva_1/chfal/virus_analyses/8054_final            # changes group of all files in listed directory
```


## 8054 Metaviral Bam

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=bwa                          # job name for listing in queue
#SBATCH --mem=50G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

echo "load any Amarel modules that script requires"
module purge # clears out any pre-existing modules
module load samtools # load any modules needed
module load bwa


ROUND_1=/projects/f_geneva_1/chfal/virus_analyses/aug_23_sequencing/round_1
ROUND_2=/projects/f_geneva_1/chfal/virus_analyses/aug_23_sequencing/round_2
# FILE=$1

echo "Bash commands for the analysis you are going to run"

echo "##################### index and align with BWA"
bwa index /projects/f_geneva_1/chfal/virus_analyses/8054_meta/8054_meta_first_scaffold.fasta

bwa mem -t 10 -5SPM /projects/f_geneva_1/chfal/virus_analyses/8054_meta/8054_meta_first_scaffold.fasta \
<(cat /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R1.fq /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R1.fq /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R1.fq) \
<(cat /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R2.fq /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R2.fq /projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R2.fq) | \
samtools sort -@10 -o /projects/f_geneva_1/chfal/virus_analyses/8054_meta/mapped_8054.bam -
samtools faidx /projects/f_geneva_1/chfal/virus_analyses/8054_meta/8054_meta_first_scaffold.fasta

echo "change user group of files created"
chgrp -R g_geneva_1 /projects/f_geneva_1/chfal/virus_analyses/8054_meta            # changes group of all files in listed directory
```


## 8158 Metaviral Bam

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --account=general
#SBATCH --job-name=bwa                          # job name for listing in queue
#SBATCH --mem=50G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=3-00:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

echo "load any Amarel modules that script requires"
module purge # clears out any pre-existing modules
module load samtools # load any modules needed
module load bwa

echo "Bash commands for the analysis you are going to run"

echo "##################### index and align with BWA"
bwa index /projects/f_geneva_1/chfal/virus_analyses/8158_meta/8158_meta_large_scaffold.fasta

bwa mem -t 10 -5SPM /projects/f_geneva_1/chfal/virus_analyses/8158_meta/8158_meta_large_scaffold.fasta \
<(cat /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R1.fq /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R1.fq) \
<(cat /projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R2.fq /projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R2.fq) | \
samtools sort -@10 -o /projects/f_geneva_1/chfal/virus_analyses/8158_meta/mapped_8158_meta.bam -
samtools faidx /projects/f_geneva_1/chfal/virus_analyses/8158_meta/8158_meta_large_scaffold.fasta

echo "change user group of files created"
chgrp -R g_geneva_1 /projects/f_geneva_1/chfal/virus_analyses/8158_meta/             # changes group of all files in listed directory
```

Now it's time to run coverage depth to calculate coverage depth on the BAM files. This is a custom script that you can just use the argument for the file which is easy to set off.

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=depth_breadth                      # job name for listing in queue
#SBATCH --mem=20G                              # memory to allocate in Mb
#SBATCH -n 10                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=02:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs

echo "load any Amarel modules that script requires"
module purge                                    # clears out any pre-existing modules
module load java
module load samtools                            # load any modules needed

echo "Bash commands for the analysis you are going to run"

echo "#...SAMTOOLS Summary stats"

samtools depth -a ${1} | awk '{c++;s+=$3}END{print s/c}' > ${1}_read_depth.txt
samtools depth -a ${1} | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > ${1}_breadth_coverage.txt
samtools flagstat ${1} > ${1}_prop_reads.txt

echo "This is a run"
echo "Now it is done"
```

We then make a table to evaluate assembly quality and decide which one we want to go with. Luckily, it seemed like including or not including the error correction step didn't really have that much of an impact on the outcome of the assembler, so I guess it doesn't really matter!

<img width="1521" height="223" alt="image" src="https://github.com/user-attachments/assets/097b4769-e11f-4a85-86bf-e27bb1cb599e" />
