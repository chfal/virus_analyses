Okay so we had to do some kind of weird awk/bash stuff to get this to work.

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

```
