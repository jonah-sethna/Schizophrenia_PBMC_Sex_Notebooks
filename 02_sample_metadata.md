---
title: "Building The Complete Sample Metadata Table"
author: "Jonah"
date: "today's date"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: show
---

# Notebook Goal

The goal  is to build a metadata table for the schizophrenia RNA-seq dataset. I will combine the original sample metadata with the validated sex assignments generated in Notebook 1.

Metadata table definition: Supplies the context needed to interpret what each sample is and how it relates to other samples.

The main inputs are:
- The counts matrix
- The original metadata
- The sex assignments

```{r}
counts <- readRDS("results/dataset/counts.rds")
samples <- read.csv("config/scz_samples_metadata.csv")
sex_assignments <- read.csv("results/sex_assignments.csv")
```


Successfully uploaded the stuff.


```{r}
head(samples)
dim(samples)
colnames(samples)
str(samples)
```

Head = First 6 Rows
Dim = Number of Rows x Columns
colnames = names all variables
str = data type of variable

Geo_accession contains a identifier used to identify each sample
Paper_sample = I don't know?
Sra_bioproject = I don't know? 
Phenotype = tells whether the sample is schizophrenia or control. 
Sample = combines the study sample name with the GEO accession to create a unique sample identifier.

```{r}
table(substr(samples$sample, 1, 3), samples$phenotype)
```
```{r}
head(colnames(counts))
head(samples$sample)
```

The phenotype labels  match the prefixes in the sample names. All 97 samples are labeled as controls, and all 84 samples are labeled as schizophrenia.

The sample column in the metadata matches the column names of the counts matrix. Sample can be used to match the metadata to the RNA-seq counts matrix.

```{r}
dim(sex_assignments)
head(sex_assignments)
colnames(sex_assignments)
```
The sex assignments table contains the sample identifier, the final validated sex assignment, and the condition. For the merge with the original metadata, I only need the sample and sex columns because phenotype information is already included in the original metadata.

```{r}
sex_clean <- sex_assignments[, c("sample", "sex")]

head(sex_clean)
```

```{r}
nrow(samples)
nrow(sex_clean)
ncol(counts)
```

```{r}
all(sex_clean$sample %in% samples$sample)

all(samples$sample %in% sex_clean$sample)
```

There are 181 samples in the original metadata, 181 samples in the validated sex assignments, and 181 sample columns in the counts matrix. All samples in the sex assignments are present in the metadata and all samples in the metadata are present in the sex assignments. 

```{r}
samples_complete <- merge(samples, sex_clean, by = "sample")
```

```{r}
head(samples_complete)
dim(samples_complete)
colnames(samples_complete)
```

The merge combined the original sample metadata with the validated sex assignments using the sample column. The merge also changed the row order.

```{r}
dim(samples_complete)
colnames(samples_complete)

nrow(samples_complete) == ncol(counts)
```

```{r}
table(samples_complete$sex, samples_complete$phenotype)
```
The completed table contains 67 female controls, 30 male controls, 46 female schizophrenia samples, and 38 male schizophrenia samples. This breakdown matches the validated sex and condition assignments from Notebook 1.

```{r}
all(samples_complete$sample == colnames(counts))
```

The order check returned FALSE, meaning that the rows of samples_complete are not in the same order as the columns of the counts matrix. Merge() reordered the rows when it combined the two tables. Need to find out how to match them

Need to fix this?

```{r}
match(colnames(counts), samples_complete$sample)
```



The match() output returned the row positions in samples_complete that correspond to each sample in the counts matrix.

```{r}
samples_complete <- samples_complete[
  match(colnames(counts), samples_complete$sample),
]
```
```{r}
all(samples_complete$sample == colnames(counts))
```

After reordering samples_complete using match(), I checked the sample order again. The comparison returned TRUE, confirming that the rows of the metadata table now exactly match the columns of the counts matrix.


How to save directory to config?


###Checkpoint

So far, I have loaded the  datasets, explored the data, validated the sample identifiers, cleaned the sex assignments, merged the tables, and checked the merged data