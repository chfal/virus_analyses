
Okay so we had our orthogroups in protein format and then we tried to run OrthoFinder in DNA sequence mode but it made totally different orthogroups, likely due to the fact that the genetic code is redundant and all that. So we basically had to make the same orthogroups from the nucleotide sequences by hand using the OrthoFinder output and the nucleotide sequences. ChatGPT was used for this because although I could have done it in R that would have required a lot of parsing so I just did this since it was faster.

First we concatenated all the VGAS nucleotide into one large fasta file.

Then we did this chatgpt code to normalize the headers so it would match with the output of OrthoFinder TSV.

```
awk '
/^>/ {
    header=$0
    sub(/^>/,"",header)

    if (match(header, /^([^ ]+ [0-9]+_[0-9]+\.\.[0-9]+)/, m))
        print ">" m[1]
    next
}
{
    print
}' all_nucleotide.fasta > normalized.fasta
```

We got the orthogroups.txt file and we turned it into a long format so we could pick through:


```
awk '
{
    og=$1
    sub(/:/,"",og)

    for(i=2;i<=NF;i++) {
        if ($i ~ /^[0-9]+$/ && (i+1)<=NF) {
            id=$i
            j=i+1

            if ($j ~ /^[0-9]+_[0-9]+\.\.[0-9]+_/) {
                split($j,a,"_")
                seq=a[1]"_"a[2]
                print og"\t"id"|"seq
                i=j
            }
        }
        else if ($i ~ /^[a-zA-Z_]+$/ && (i+1)<=NF) {
            id=$i
            j=i+1

            if ($j ~ /^[0-9]+_[0-9]+\.\.[0-9]+_/) {
                split($j,a,"_")
                seq=a[1]"_"a[2]
                print og"\t"id"|"seq
                i=j
            }
        }
    }
}' orthogroups.txt > orthogroup_sequences.tsv
```

We're reformatting to a kind of key value pair lookup thing:
```
awk -F'\t' '{print $2 "\t" $1}' orthogroup_sequences.tsv > seq_to_og.txt
```

Now we are making the orthogroups and they look great!
```
awk '
BEGIN { FS="\t" }

# Read: sequence_ID -> orthogroup
FNR==NR {
    og[$1] = $2
    next
}

# FASTA header
/^>/ {
    header = substr($0,2)

    if (header in og) {
        current_og = og[header]
        outfile = "orthogroup_fastas/" current_og ".fasta"

        print $0 >> outfile
        keep = 1
    }
    else {
        keep = 0
    }

    next
}

# FASTA sequence
keep {
    print $0 >> outfile
}
' seq_to_og.txt normalized.fasta
```

So now we have all the fastas together grouped but because they are not aligned we need to align them to each other:

```
#!/bin/bash
#SBATCH --partition=p_ccib_1
#SBATCH --account=general
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=masce_10
#SBATCH --mem=20G
#SBATCH -n 15
#SBATCH -N 1
#SBATCH --time=7-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL,END,REQUEUE



module purge
module load java

for fasta in *.fasta; do
    base="${fasta%.fasta}"

    java -jar macse_v2.07.jar \
        -prog alignSequences \
        -seq "$fasta" \
        -out_NT "${base}_nt_aligned.fasta" \
        -out_AA "${base}_aa_aligned.fasta"
done
```

Now lets play this game where we remove stop codons from the alignments

```
 for file in *.fasta; do     base=$(basename "$file" .fasta)
    ./hyphy CleanStopCodons.bf <<EOF$file
${base}_cleaned.fasta
EOF
 done

```

Then we will do iqtree. we are running this without any premonition of what is in the foreground or the background.


```
#!/bin/bash
#SBATCH --partition=p_geneva_1
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=iqtree
#SBATCH --mem=10G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


module purge
eval "$(conda shell.bash hook)"
conda activate iqtree



for file in *nt_cleaned.fasta; do
    base="${file%.fasta}"
    iqtree -s "$file" -m MFP -bb 1000 -alrt 1000 -pre "$base"
done
```

finally: absrel

```
#!/bin/bash
#SBATCH --partition=p_geneva_1
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=iqtree
#SBATCH --mem=10G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL

module purge
eval "$(conda shell.bash hook)"
conda activate hyphy


mkdir -p absrel_results

for file in *_nt_cleaned.fasta; do
    base="${file%.fasta}"

    hyphy absrel \
        --alignment "$file" \
        --tree "${base}.treefile" \
        --output "absrel_results/${base}.json" \
        | tee -a "absrel_results/${base}.log"
done
```

# but wait!
# realigning (chakras) once again

I forgot with cristatellus we actually redid the alignments using another version of macse so I just did that really quickly.

```
#!/bin/bash
#SBATCH --partition=p_geneva_1                    # which partition to run the job, options are in the Amarel guide
#SBATCH --exclude=gpuc001,gpuc002               # exclude CCIB GPUs
#SBATCH --job-name=macse_alfix                      # job name for listing in queue
#SBATCH --mem=10G                              # memory to allocate in Mb
#SBATCH -n 1                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --cpus-per-task=24
#SBATCH --time=3-00:00:00                         # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                               # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=FAIL,END,REQUEUE      # email for the following reasons


module purge
module load apptainer/1.2.5-sg1509

# Had to run it like this otheriwse it wouldnt find the files needed for the pipeline to work

for fasta in *.fasta; do
    base="${fasta%.fasta}"

    echo "Processing $fasta..."

    apptainer run \
        --bind /projects/f_geneva_1/chfal:/projects/f_geneva_1/chfal \
        --pwd /projects/f_geneva_1/chfal/virus_analyses/genes/orthogroup_nucleotides \
        MACSE_ALFIX_v01.sif \
        --out_dir . \
        --out_file_prefix "${base}_aligned" \
        --in_seq_file "$fasta"
done


```

And then we do IQTREE / ABSREL

iqtree: 

```
#!/bin/bash
#SBATCH --partition=p_geneva_1
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=iqtree
#SBATCH --mem=10G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL


module purge
eval "$(conda shell.bash hook)"
conda activate iqtree



for file in *aln; do
    base="${file%.aln}"
    iqtree -s "$file" -m MFP -bb 1000 -alrt 1000 -pre "$base"
done

```


Absrel: 

```
#!/bin/bash
#SBATCH --partition=p_geneva_1
#SBATCH --exclude=gpuc001,gpuc002
#SBATCH --job-name=iqtree
#SBATCH --mem=10G
#SBATCH -n 10
#SBATCH -N 1
#SBATCH --time=3-00:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=FAIL

module purge
eval "$(conda shell.bash hook)"
conda activate hyphy


mkdir -p absrel_results

for file in *NT_cleaned.aln; do
    base="${file%.aln}"

    hyphy absrel \
        --alignment "$file" \
        --tree "${base}.treefile" \
        --output "absrel_results/${base}.json" \
        | tee -a "absrel_results/${base}.log"
done

```

We ran OpenRDP on only the alignments that passed macse alfix because the plain alignments from macse were not working as they had lots of question marks / kind of were bad alignments. This makes a bunch of .csv files I truly have no clue how to parse.

```
for file in *_cleaned.aln; do
    base="${file%_cleaned.aln}"
    openrdp "$file" -m rdp -o "${base}_RDP.csv"
done

```
