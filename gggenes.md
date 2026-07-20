# Gggenes

Gggenes will allow us to visualize the identities of genes and their locations on the chromosome.

However we need to prepare our input files for gggenes.

We had to go back to the *annotation_out.out file from vgas, which is a table that provides a bunch of information about the predicted genes for each molecule.


```
# this code removes the first 12 lines of the file which are not table-=shaped and just have random information in them

 for f in *; do     sed -i '1,12d' "$f"; done

# this renames the files and removes the _annotation part of the file name which helps in the later step of making the table

for f in *_annotation.out; do
    mv "$f" "${f/_annotation/}"
done


# this takes the file name and appends it to column 1

for f in *.out; do
    prefix="${f%.out}"

    awk -v p="$prefix" '
    BEGIN{OFS="\t"}

    {
        $1 = p "_" $1
        print
    }
    ' "$f" > tmp && mv tmp "$f"
done
```


i then put this all into excel and did a little bit of excel magic to get this:

```
molecule	gene	start	end	strand	orientation
8054	8054_1	79	855	+	
8054	8054_2	902	1042	-	
8054	8054_3	1309	2739	+	
8054	8054_4	2767	4077	+	
8054	8054_5	4140	5492	+	
8054	8054_6	5551	6933	+	
8054	8054_7	6965	8368	+	
8054	8054_8	8433	9068	-	
8054	8054_9	9019	9450	-	
8054	8054_10	9456	9635	+	
```

However we want to have blast identities too so I had to blast them and I mostly did it in excel with just cuttiing and pasting but i got a column that looked like this:

<img width="258" height="609" alt="image" src="https://github.com/user-attachments/assets/0b5492ef-5439-4ab9-bf22-67388202a104" />


and then i left joined them together in R with left join (see file gggenes.R)
