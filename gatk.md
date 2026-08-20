# Running GATK on viruses:

Run BWA individually on the different BAM files:

```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --account=general
#SBATCH --job-name=bwa
#SBATCH --mem=50G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

echo "load any Amarel modules that script requires"
module purge
module load samtools
module load bwa

ROUND_1=/projects/f_geneva_1/chfal/virus_analyses/aug_23_sequencing/round_1
ROUND_2=/projects/f_geneva_1/chfal/virus_analyses/aug_23_sequencing/round_2

REF=/projects/f_geneva_1/chfal/virus_analyses/gatk/8054.fasta
OUTDIR=/projects/f_geneva_1/chfal/virus_analyses/gatk

echo "##################### index reference"
bwa index $REF
samtools faidx $REF

echo "##################### map sequencing run 1"

bwa mem -t 10 $REF \
/projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R1.fq \
/projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8054_S1_L001_filtered.R2.fq \
| samtools sort -@10 -o $OUTDIR/mapped_8054_run1.bam -

samtools index $OUTDIR/mapped_8054_run1.bam


echo "##################### map sequencing run 2"

bwa mem -t 10 $REF \
/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R1.fq \
/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8054_S1_L001_filtered.R2.fq \
| samtools sort -@10 -o $OUTDIR/mapped_8054_run2.bam -

samtools index $OUTDIR/mapped_8054_run2.bam


echo "##################### map sequencing run 3"

bwa mem -t 10 $REF \
/projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R1.fq \
/projects/f_geneva_1/chfal/virus_analyses/8054_2026/8054_S1_L001_filtered.R2.fq \
| samtools sort -@10 -o $OUTDIR/mapped_8054_run3.bam -

samtools index $OUTDIR/mapped_8054_run3.bam


echo "##################### change user group of files created"

chgrp -R g_geneva_1 $OUTDIR

echo "##################### DONE"
echo "Created:"
echo "$OUTDIR/mapped_8054_run1.bam"
echo "$OUTDIR/mapped_8054_run2.bam"
echo "$OUTDIR/mapped_8054_run3.bam"

```


```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --account=general
#SBATCH --job-name=bwa
#SBATCH --mem=50G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

echo "load any Amarel modules that script requires"
module purge
module load samtools
module load bwa

REF=/projects/f_geneva_1/chfal/virus_analyses/gatk/8158.fasta
OUTDIR=/projects/f_geneva_1/chfal/virus_analyses/gatk

echo "##################### index reference"
bwa index $REF
samtools faidx $REF

echo "##################### map sequencing run 1"

bwa mem -t 10 $REF \
/projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R1.fq \
/projects/f_geneva_1/data/virus_sequencing/20230130_FS10002629_1_BRR99421-1634/Alignment_1/20230131_080321/Fastq/8158_S2_L001_filtered.R2.fq \
| samtools sort -@10 -o $OUTDIR/mapped_8158_run1.bam -

samtools index $OUTDIR/mapped_8158_run1.bam


echo "##################### map sequencing run 2"

bwa mem -t 10 $REF \
/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R1_v2.fq \
/projects/f_geneva_1/data/virus_sequencing/20230215_FS10002629_2_BRR99421-2334/Alignment_1/20230216_031252/Fastq/8158_S2_L001_filtered.R2_v2.fq \
| samtools sort -@10 -o $OUTDIR/mapped_8158_run2.bam -

samtools index $OUTDIR/mapped_8158_run2.bam


echo "##################### change user group of files created"

chgrp -R g_geneva_1 $OUTDIR

echo "##################### DONE"
echo "Created:"
echo "$OUTDIR/mapped_8158_run1.bam"
echo "$OUTDIR/mapped_8158_run2.bam"

```


Add replace read groups, still individually:

```

#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002,halc068
#SBATCH --job-name=addorreplace
#SBATCH --mem=50G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL

echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load gatk
module load samtools

echo "load variables"

SAMPLE=8054
OUTDIR="/projects/f_geneva_1/chfal/virus_analyses/gatk"
BAM_DIR="/projects/f_geneva_1/chfal/virus_analyses/gatk"

echo "##################### Add read group to run 1"

BAM=${BAM_DIR}/mapped_8054_run1.bam
PU=BRR99421-1634:1:1101

gatk AddOrReplaceReadGroups \
-I ${BAM} \
-O ${OUTDIR}/mapped_8054_run1.addRG.bam \
-LB library1 \
-PL illumina \
-PU ${PU} \
-SM ${SAMPLE}

samtools index ${OUTDIR}/mapped_8054_run1.addRG.bam


echo "##################### Add read group to run 2"

BAM=${BAM_DIR}/mapped_8054_run2.bam
PU=BRR99421-2334:1:1101

gatk AddOrReplaceReadGroups \
-I ${BAM} \
-O ${OUTDIR}/mapped_8054_run2.addRG.bam \
-LB library1 \
-PL illumina \
-PU ${PU} \
-SM ${SAMPLE}

samtools index ${OUTDIR}/mapped_8054_run2.addRG.bam


echo "##################### Add read group to run 3"

BAM=${BAM_DIR}/mapped_8054_run3.bam
PU=BWB90213-2306:1:1101

gatk AddOrReplaceReadGroups \
-I ${BAM} \
-O ${OUTDIR}/mapped_8054_run3.addRG.bam \
-LB library1 \
-PL illumina \
-PU ${PU} \
-SM ${SAMPLE}

samtools index ${OUTDIR}/mapped_8054_run3.addRG.bam


echo "##################### change user group"

chgrp -R g_geneva_1 ${OUTDIR}

echo "##################### DONE"

echo "Created:"
echo "${OUTDIR}/mapped_8054_run1.addRG.bam"
echo "${OUTDIR}/mapped_8054_run2.addRG.bam"
echo "${OUTDIR}/mapped_8054_run3.addRG.bam"
```


```
#!/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=gpuc001,gpuc002,halc068
#SBATCH --job-name=addorreplace
#SBATCH --mem=50G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL

echo "load modules"
module purge
module use /projects/community/modulefiles/
module load java
module load gatk
module load samtools

echo "load variables"

SAMPLE=8158
OUTDIR="/projects/f_geneva_1/chfal/virus_analyses/gatk"
BAM_DIR="/projects/f_geneva_1/chfal/virus_analyses/gatk"

echo "##################### Add read group to run 1"

BAM=${BAM_DIR}/mapped_8158_run1.bam
PU=BRR99421-1634:1:1101

gatk AddOrReplaceReadGroups \
-I ${BAM} \
-O ${OUTDIR}/mapped_8158_run1.addRG.bam \
-LB library1 \
-PL illumina \
-PU ${PU} \
-SM ${SAMPLE}

echo "index run 1"
samtools index ${OUTDIR}/mapped_8158_run1.addRG.bam


echo "##################### Add read group to run 2"

BAM=${BAM_DIR}/mapped_8158_run2.bam
PU=BRR99421-2334:1:1101

gatk AddOrReplaceReadGroups \
-I ${BAM} \
-O ${OUTDIR}/mapped_8158_run2.addRG.bam \
-LB library1 \
-PL illumina \
-PU ${PU} \
-SM ${SAMPLE}

echo "index run 2"
samtools index ${OUTDIR}/mapped_8158_run2.addRG.bam


echo "##################### DONE"

# read group going with reads from 1-30
# BRR99421-1634:1:1101

# read group from 2-15
# BRR99421-2334:1:1101

```

Merge for both 8054/8158 separately:


```
gatk MergeSamFiles -I mapped_8054_run1.addRG.bam -I mapped_8054_run2.addRG.bam -I mapped_8054_run3.addRG.bam -O merged_8054.addRG.bam

```


```
gatk MergeSamFiles -I mapped_8158_run1.addRG.bam -I mapped_8158_run2.addRG.bam -O merged_8158.addRG.bam

```

Mark duplicates:

```
gatk MarkDuplicates \
-I merged_8054.addRG.bam \
-O merged_8054.addRG.marked.bam \
-M 8054_merged_metrics.txt \
--REMOVE_DUPLICATES false --ASSUME_SORTED true --CREATE_INDEX true
```

```
gatk MarkDuplicates \
-I merged_8158.addRG.bam \
-O merged_8158.addRG.marked.bam \
-M 8054_merged_metrics.txt \
--REMOVE_DUPLICATES false --ASSUME_SORTED true --CREATE_INDEX true

```


Run haplotype caller with ploidy of 1 for virus:

```
gatk HaplotypeCaller \
--native-pair-hmm-threads 2 \
-I merged_8054.addRG.marked.bam \
-O 8054.g.vcf.gz \
-R 8054.fasta \
-ERC BP_RESOLUTION \
--output-mode EMIT_ALL_CONFIDENT_SITES \
--max-reads-per-alignment-start 0 \
--sample-ploidy 1 \
-RF NotDuplicateReadFilter
```

```

gatk HaplotypeCaller \
--native-pair-hmm-threads 2 \
-I merged_8158.addRG.marked.bam \
-O 8158.g.vcf.gz \
-R 8158.fasta \
-ERC BP_RESOLUTION \
--output-mode EMIT_ALL_CONFIDENT_SITES \
--max-reads-per-alignment-start 0 \
--sample-ploidy 1 \
-RF NotDuplicateReadFilter
```
