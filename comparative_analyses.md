# Expanding Analyses to Other Atadenoviruses

 We downloaded a bunch of ATDV sequences from GenBank (check the excel spreadsheet supplement / supplemental table to learn which ones). We had a total of 35 downloaded sequences.

Simple script to rename all *.fna downloaded from GenBank to .fasta

```
for f in *.fna; do
    mv -- "$f" "${f%.fna}.fasta"
done
```

Now we should annotate these with VGAS using the regular pipeline script to do it. 

```
!#/bin/bash
#SBATCH --partition=cmain
#SBATCH --exclude=halc068
#SBATCH --job-name=vgas
#SBATCH --mem=10G
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --time=5:00:00
#SBATCH --requeue
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu
#SBATCH --mail-type=END,FAIL

module purge
eval "$(conda shell.bash hook)"
conda activate vgas

for fasta in *.fasta; do

    base=$(basename "$fasta" .fasta)

    echo "Running $base"

    vgas "$fasta" "${base}_annotation.fasta" -p

done
```

We also have to do the thing where we cut the first four lines of the fasta and remove any white spaces in each protein fasta.

```

sed -i '1,4d' protein*.fasta
sed -i 's/:/_/g' *.fasta
```

Next we'll reformat these fasta files

```
for f in protein_*_annotation.fasta; do mv "$f" "$(echo "$f" | sed 's/^protein_//; s/_annotation//')"; done
```

```
for f in *.fasta; do
    name="${f%.fasta}"
    sed -i "s/^>Potential protein/>${name}/" "$f"
done
```


Then we can run OrthoFinder

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --account=general
#SBATCH --constraint=oarc
#SBATCH --exclude=halc068               # exclude CCIB GPUs
#SBATCH --job-name=orthofinder                        # job name for listing in queue
#SBATCH --mem=10G                               # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=5:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                                # restart and paused or superseeded jobs
#SBATCH --mail-user=chf29@scarletmail.rutgers.edu           # email address to send status updates
#SBATCH --mail-type=END                     # email for the following reasons

echo "Load conda needed for orthofinder"

module purge
eval "$(conda shell.bash hook)"
conda activate orthofinder

cd /projects/f_geneva_1/chfal/virus_analyses/test_orthofinder

echo "Create variables for Orthofinder"

ulimit -n 2400

orthofinder -f /projects/f_geneva_1/chfal/virus_analyses/test_orthofinder -M msa                  # Run full OrthoFinder analysis on FASTA format proteomes in specfied directory

# orthofinder [options] -f <dir1> -b <dir2>     # Add new species in to a previous run and run new analysis
```


Now we have our orthogroup table. We have almost all of it done in genecounts.tsv except for the singletons (which is annoying.)

here we are making a two-column file of the orthogroup ID and then the species it came from.

```

 grep ">" *.fa > names.txt

cut -f 1 -d "." names.txt > orthogroup_names.txt

cut -f 2 -d ">" names.txt | cut -f 1 -d " " > species.txt

paste orthogroup_names.txt species.txt > names_species.txt


```


now we are using awk to append these onto the bottom of the orthogroup gene count file.

```
BEGIN{FS=OFS="\t"}

NR==FNR{
    og[NR]=$1
    sp[NR]=$2
    next
}

FNR==1{
    for(i=2;i<=NF;i++){
        if($i=="Total") continue
        col[$i]=i
    }
    n=NF
    print $0
    next
}

{
    print $0
    next
}

END{
    for(i=1;i<=length(og);i++){

        printf og[i]

        total=0

        for(j=2;j<n;j++){
            if(j==col[sp[i]]){
                printf OFS 1
                total=1
            } else {
                printf OFS 0
            }
        }

        printf OFS total "\n"
    }
}
' names_species.tsv Orthogroups.GeneCount.tsv > Orthogroups.GeneCount.All.tsv

```


then i copied that file into a new file with transpose (i wanted columns to become rows and rows to become columns so i could add the columns of the clades on the Y axis.


okay so now we need to remake the tree again. i decided to use these orthogroups because we wanted to have orthogroups that were in more than 40% of our species which this is the cutoff (the remaining orthogroups have 9 or below sequences, there is a very deep dropoff here):

```
OG0000000.fa:39
OG0000001.fa:37
OG0000002.fa:37
OG0000003.fa:37
OG0000004.fa:37
OG0000005.fa:37
OG0000006.fa:37
OG0000007.fa:36
OG0000008.fa:36
OG0000009.fa:36
OG0000010.fa:36
OG0000011.fa:35
OG0000012.fa:35
OG0000013.fa:30
OG0000014.fa:29
OG0000015.fa:23
OG0000016.fa:23
OG0000017.fa:23
OG0000018.fa:20
OG0000019.fa:18
OG0000020.fa:18
OG0000021.fa:16
OG0000022.fa:16

```

It turns out 8054 and 8158 are present in all 23 of these orthogroups.

Below is the list of orthogroups that 8054 and 8158 are present in.

```
OG0000000.fa:>8054_16_13931..14977__348_aa
OG0000001.fa:>8054_23_23025..23663__212_aa
OG0000002.fa:>8054_26_24406..25788__460_aa
OG0000003.fa:>8054_27_25822..27651__609_aa
OG0000004.fa:>8054_28_27608..28615__335_aa
OG0000005.fa:>8054_31_28996..30429__477_aa
OG0000006.fa:>8054_32_30462..33665__1067_aa
OG0000007.fa:>8054_19_16204..18348__714_aa
OG0000008.fa:>8054_20_18388..19659__423_aa
OG0000009.fa:>8054_21_19681..20286__201_aa
OG0000010.fa:>8054_22_20283..23006__907_aa
OG0000011.fa:>8054_18_15176..15844__222_aa
OG0000012.fa:>8054_33_33371..34723__450_aa
OG0000013.fa:>8054_15_13218..13898__226_aa
OG0000014.fa:>8054_25_24000..24365__121_aa
OG0000015.fa:>8054_12_11133..12029__298_aa
OG0000016.fa:>8054_13_12150..12599__149_aa
OG0000017.fa:>8054_14_12590..13270__226_aa
OG0000018.fa:>8054_34_34753..35865__370_aa
OG0000019.fa:>8054_24_23745..23996__83_aa
OG0000020.fa:>8054_36_36311..37285__324_aa
OG0000021.fa:>8054_17_14995..15165__56_aa
OG0000022.fa:>8054_35_35911..36312__133_aa
OG0000024.fa:>8054_4_2767..4077__436_aa
OG0000024.fa:>8054_6_5551..6933__460_aa
OG0000024.fa:>8054_7_6965..8368__467_aa
OG0000026.fa:>8054_3_1309..2739__476_aa
OG0000026.fa:>8054_5_4140..5492__450_aa
OG0000034.fa:>8054_11_9767..11038__423_aa
OG0000054.fa:>8054_1_79..855__258_aa
OG0000055.fa:>8054_2_902..1042__46_aa
OG0000056.fa:>8054_8_8433..9068__211_aa
OG0000057.fa:>8054_9_9019..9450__143_aa
OG0000058.fa:>8054_10_9456..9635__59_aa
OG0000059.fa:>8054_29_28605..28763__52_aa
OG0000060.fa:>8054_30_28773..28916__47_aa


(vgas) [chf29@amarel1 Orthogroup_Sequences]$ grep "8158" *.fa
OG0000000.fa:>8158_16_13940..14986__348_aa
OG0000001.fa:>8158_23_23034..23672__212_aa
OG0000002.fa:>8158_26_24415..25797__460_aa
OG0000003.fa:>8158_27_25831..27660__609_aa
OG0000004.fa:>8158_28_27617..28624__335_aa
OG0000005.fa:>8158_31_29005..30438__477_aa
OG0000006.fa:>8158_32_30471..33674__1067_aa
OG0000007.fa:>8158_19_16213..18357__714_aa
OG0000008.fa:>8158_20_18397..19668__423_aa
OG0000009.fa:>8158_21_19690..20295__201_aa
OG0000010.fa:>8158_22_20292..23015__907_aa
OG0000011.fa:>8158_18_15185..15853__222_aa
OG0000012.fa:>8158_33_33380..34732__450_aa
OG0000013.fa:>8158_15_13227..13907__226_aa
OG0000014.fa:>8158_25_24009..24374__121_aa
OG0000015.fa:>8158_12_11142..12038__298_aa
OG0000016.fa:>8158_13_12159..12608__149_aa
OG0000017.fa:>8158_14_12599..13279__226_aa
OG0000018.fa:>8158_34_34762..35874__370_aa
OG0000019.fa:>8158_24_23754..24005__83_aa
OG0000020.fa:>8158_36_36320..37294__324_aa
OG0000021.fa:>8158_17_15004..15174__56_aa
OG0000022.fa:>8158_35_35920..36321__133_aa
OG0000024.fa:>8158_4_2776..4086__436_aa
OG0000024.fa:>8158_6_5560..6942__460_aa
OG0000024.fa:>8158_7_6974..8377__467_aa
OG0000026.fa:>8158_3_1318..2748__476_aa
OG0000026.fa:>8158_5_4149..5501__450_aa
OG0000034.fa:>8158_11_9776..11047__423_aa
OG0000054.fa:>8158_1_88..864__258_aa
OG0000055.fa:>8158_2_911..1051__46_aa
OG0000056.fa:>8158_8_8442..9077__211_aa
OG0000057.fa:>8158_9_9028..9459__143_aa
OG0000058.fa:>8158_10_9465..9644__59_aa
OG0000059.fa:>8158_29_28614..28772__52_aa
OG0000060.fa:>8158_30_28782..28925__47_aa
```



As decided before we wanted to use the first 23 orthogroups which included the majority of the species in OrthoFinder (OG00000-OG000023).

We then had to run IQTree on all of those multiple sequence alignments, which was done like this in one single lazy-ass for loop:

```
#!/bin/bash
#SBATCH --partition=cmain                    # which partition to run the job, options are in the Amarel guide
#SBATCH --account=general
#SBATCH --exclude=halc068               # exclude CCIB GPUs
#SBATCH --job-name=iqtree                        # job name for listing in queue
#SBATCH --mem=10G                               # memory to allocate in Mb
#SBATCH -n 20                                   # number of cores to use
#SBATCH -N 1                                    # number of nodes the cores should be on, 1 means all cores on same node
#SBATCH --time=5:00:00                       # maximum run time days-hours:minutes:seconds
#SBATCH --requeue                                # restart and paused or superseeded jobs

echo "Load conda needed for orthofinder"

module purge
eval "$(conda shell.bash hook)"
conda activate iqtree

for FILE in *.fa
do
    echo "Running IQ-TREE on ${FILE}"
        iqtree -m TEST -s ${FILE} -redo

done

```

We got then a bunch of tree files out, which we concatenated into a single file called input tree:

```
cat *treefile > input_tree.tre
```


Then it was time to use the new version of Astral, which is actually called Aster now:


https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral-pro3.md

I could have run this on the cluster but was running into some compatibility issues so what I did was I installed the Windows application packages and that worked great. I also had to trim the tree labels because Aster will accept multiple phylogenetic trees with both duplicate values and incomplete taxa (since not all taxa are represented in this). However there were lots of fragmented names left over from running IQTREE, so I stripped the out using this perl command (chatgpt helped here).

```
perl -pe 's/([A-Za-z0-9]+(?:_[A-Za-z0-9]+)*)_\1(?:_\d+_\d+\.\.\d+_)?/$1/g;' input.tree > renamed.tree
```

Once I had that, I had a renamed.tree file that I clicked and ran through the Astral-Pro GUI and it produced this lovely tree:

<img width="1037" height="953" alt="image" src="https://github.com/user-attachments/assets/f959ae7e-a6f5-4879-ad43-61262ace8923" />


Which is going to be the input of the next graph I will try to make!


# Running IQTREE again:

We want to find out which genes are under which signatures of selection

We renamed all the Adenovirus files in the alignment from OG000000-OG00000023.fa to just be the species names and not the protein IDs and we concatenated them into a concatenated nexus and phylip file (attached)

We downloaded IQTree and ran this in an interactive job

```
/iqtree3 -s sequences.phy -p nexus.nex -m TEST -bb 1000 -redo
```
