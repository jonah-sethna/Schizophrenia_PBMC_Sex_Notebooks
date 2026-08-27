---
title: "Exploring a Bulk RNA-Seq Counts Matrix"
author: "Jonah"
date: "August 19, 2026"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: show
  pdf_document:
    toc: true
---

## Part 1: Loading and getting to know the data



I am first loading the counts matrix and the features table into R



```{r}
counts <- readRDS("results/dataset/counts.rds")
features <- readRDS("results/dataset/features.rds")
```
---



Both the counts matrix and the features table loaded successfully without an error.



### Checking the object types



I am checking what type of R object was created from each file.



```{r}
class(counts)
class(features)
```



The counts object is a matrix and array. Need to ask Helena what that means. GPT says that it just means it is 2-D structure with rows and columns.The features table comes back data.frame. This may mean it is a it describes the features of the matrix rows.


### Checking the size of the counts matrix



I am checking the dimesniosn of the counts matrix to determine how many feautures and samples it contains.




```{r}
dim(counts)
```



Seems like 89118 rows and 181 columns.
The rows are the Genes/HERVS and the columns are individual samples



### Previewing the counts matrix



I guess this displays the first 5 rows and columns to see the matrix and how the data looks


```{r}
counts[1:5, 1:5]
```



This is pretty cool. On the top we have the indivudal people and the left is probably the HERVS or genes or whatever.
The first sections of the matrix contains mostly zero or low counts. One feature had counts between 4-7 which is the highest compared to other samples.
Most features had zero counts in their samples.
This shows that features may not be detected in many samples.



### Examining the feature names



I am examining the first and last 20 row names to look for patterns in how they are labeled



```{r}
head(rownames(counts), 20)
tail(rownames(counts), 20)
```



Row names with HERV before it seem to be HERVS. I am not familiar with what LINES and CG look like. 

> features$gene.type[953]
[1] "CG"
> features$gene.type[which(features$gene.name == "DNAJC8")]
[1] "CG"
> features$gene.type[which(features$gene.name == "DNAJC8"), ]
Error in features$gene.type[which(features$gene.name == "DNAJC8"), ] : 
  incorrect number of dimensions

> features[which(features$gene.name == "DNAJC8"), ]
    gene.name gene.type
953    DNAJC8        CG
> features[953, ]
    gene.name gene.type
953    DNAJC8        CG



```{r}
head(rownames(counts))
tail(rownames(counts))
```


```{r}
#rownames(counts)
```





### Examining the samples names

```{r}
colnames(counts)
```




SCZ = Schizophrenia
Con = Healthy
GSM = ???




```{r}
ncol(counts)

grepl("SCZ", colnames(counts))

sum(grepl("SCZ", colnames(counts)))
sum(grepl("CON", colnames(counts)))
```
The matrix contains 181 samples. 84 schizophrenia samples and 97 control samples.




```{r}
head(features)
tail(features)
dim(features)
table(features$gene.type)
```
```{r}
dim(features[which(features$gene.name == "DNAJC8"), ])
```


The features table classifies 60,605 features as CG, 13,545 as LINE elements, and 14,968 as LTR/HERV elements.




### Comparing the features table with the counts matrix

```{r}
nrow(features) == nrow(counts)
```




Both the features table and the counts matrix contain 89,118 rows.



## Part 2: Playing with the matrix

### Accessing a specific count

I am retrieving the RNA-seq count for the WASH7P feature to practice selecting one matrix.

```{r}
counts["WASH7P", "SCZ47_GSM7796488"]
```

The WASH7P count in the sample was 4. 

### Examining WASH7P across all samples



I am getting the WASH7P row to examine its expression across all 181 samples. Then make into bar plot.

```{r}
wash7p_counts <- counts["WASH7P", ]

wash7p_counts
mean(wash7p_counts)
max(wash7p_counts)

barplot(
  wash7p_counts,
  las = 2,
  cex.names = 0.3,
  main = "WASH7P expression across all samples",
  ylab = "Raw counts"
)
```
```{r}
mean(wash7p_counts)
```

```{r}
max(wash7p_counts)
```

WASH7P expression varied. Most samples had  low raw counts. Some had much higher expression though The mean count was 7.088, and the maximum count was 30.

### Examining DDX11L1 across all samples

Helena wants me to sample another gene and write what I see. I will compare DDX11L1 expression across all samples and its distribution with WASH7P.

```{r}
ddx11l1_counts <- counts["DDX11L1", ]

mean(ddx11l1_counts)
max(ddx11l1_counts)

barplot(
  ddx11l1_counts,
  las = 2,
  cex.names = 0.3,
  main = "DDX11L1 expression across all samples",
  ylab = "Raw counts"
)
```
```{r}
c(
  mean = mean(ddx11l1_counts),
  maximum = max(ddx11l1_counts)
)
```
DDX11L1 had a mean raw count of  0.60 and a maximum count of 7. Its expression was much lower than WASH7P.

### Examining a single sample

I am examining all 89,118 features in one sample.

```{r}
scz47 <- counts[, "SCZ47_GSM7796488"]

zero_features <- sum(scz47 == 0)
zero_fraction <- sum(scz47 == 0) / length(scz47)

c(
  zero_features = zero_features,
  zero_fraction = zero_fraction
)
```
```{r}
scz47 <- counts[, "SCZ47_GSM7796488"]

zero_features <- sum(scz47 == 0)
zero_fraction <- sum(scz47 == 0) / length(scz47)

c(
  zero_features = zero_features,
  zero_fraction = zero_fraction
)
```
```{r}
hist(
  scz47,
  breaks = 100,
  main = "Distribution of feature counts in SCZ47",
  xlab = "Raw counts"
)

hist(
  scz47[scz47 > 0 & scz47 < 1000],
  breaks = 100,
  main = "Non-zero counts in SCZ47 (< 1000)",
  xlab = "Raw counts"
)
```

In the sample 58,294 genes had zero counts which is 65.41% of all features. This indicates that most features in the matrix were not detected in this PBMC sample. This is reasonable because not every gene or transposable element is active.

```{r}
hist(
  scz47,
  breaks = 100,
  main = "Distribution of feature counts in SCZ47",
  xlab = "Raw counts"
)

hist(
  scz47[scz47 > 0 & scz47 < 1000],
  breaks = 100,
  main = "Non-zero counts in SCZ47 (< 1000)",
  xlab = "Raw counts"
)
```

Most features had zero or low raw counts, while a small number had extremely high counts reaching 120,000. 

```{r}
lib_sizes <- colSums(counts)
summary(lib_sizes)
```
```{r}
barplot(
  lib_sizes,
  las = 2,
  cex.names = 0.3,
  main = "Library size per sample",
  ylab = "Total counts"
)
```

Library Defintion - total number of RNA-seq reads counted in one sample—the sum of all values in that sample’s column.
The sample library sizes ranged from 10000000 to 22000000 total counts. Most samples had broadly similar library sizes, although several samples had noticeably higher totals.

### Calculating total counts for each feature

I am summing each feature's counts across all 181 samples. I will identify the 20 features with the highest total expression and determine how many features have zero counts in every sample.


```{r}
gene_totals <- rowSums(counts)

top_20_features <- tail(sort(gene_totals), 20)
top_20_features

sum(gene_totals == 0)
```
The most highly expressed feature across all samples was B2M. A total of 9,532 features had zero counts across every sample. These completely undetected features provide no information for comparing samples and should be removed.


## Checkpoint: What I Have Learned So Far

The dataset contains raw RNA-seq counts from PBMC samples. The rows represent genes, LINE elements, and HERV elements. The columns represent individual samples.

The features include: 60,605 canonical genes, 13,545 LINE elements, 14,968 LTR/HERV elements**

For WASH7P:

Count in sample SCZ47: 4
Mean across all samples: 7.09
Maximum: 30

For DDX11L1:

Mean across all samples: 0.60
Maximum: 7

WASH7P generally had higher raw counts than DDX11L1.

### The SCZ47 sample

In SCZ47: 58,294 features had zero counts.
This was 65.41% of all features.

The histograms showed that most features had zero or low counts, while a small number had extremely high counts.

Library size is the total number of RNA-seq counts in one sample.

* Minimum: ~ 10,000,000
* Median: ~ 14,000,000
* Mean: ~ 14,000,000
* Maximum: ~ 22,000,000

Most samples had similar library sizes.

### Total expression across all samples

`B2M` was the most highly expressed feature, with 24,090,715 total counts.

Other highly expressed features included genes such as HLA genes.

A total of 9,532 features, or approximately 10.7%, had zero counts in every sample.

### Main conclusion

We must filter out the stuff that is useless to me now.







### Free exploration

Helena wants me to do a bunch of questions. I need to now investigate how many features had high or extremely low total counts across all 181 samples.

```{r}
sum(gene_totals > 1000)
sum(gene_totals > 100)
sum(gene_totals < 10)
```
```{r}
head(sort(scz47, decreasing = TRUE), 10)
```

Of the 89,118 features, 20,325 had more than 1,000 total counts, 35,165 had more than 100 total counts, and 23,470 had fewer than 10 total counts. 

B2M was the most highly expressed feature in SCZ47, with 121,356 counts. 


```{r}
#counts["NONEXISTENT_GENE", ]
```


## Part 3: Filtering

I am examining the fraction of features with zero or very low total counts. Features with consistently low counts may provide little useful information and can add noise to later analyses. Need to filter the un-needed crap.


I am first calculating what fraction of features have total counts below different thresholds.

```{r}
gene_totals <- rowSums(counts)

mean(gene_totals == 0)
mean(gene_totals < 10)
mean(gene_totals < 50)
```

10% of features had zero total counts, 26% had fewer than 10 counts, and 50% had fewer than 50 counts across all 181 samples. 

Now I need to figure out a filtering rule:

### Creating the filtering rule

I will keep features that have at least five counts in at least ten samples because thats what Helena said.

```{r}
# TRUE means the count is at least 5
(counts >= 5)[1:5, 1:5]

# Count how many samples meet this requirement for each feature
genes_above_5 <- rowSums(counts >= 5)
head(genes_above_5, 20)

# Keep features that meet the requirement in at least 10 samples
keep <- genes_above_5 >= 10

c(
  features_kept = sum(keep),
  features_removed = sum(!keep),
  fraction_kept = mean(keep)
)
```

The filter kept 29831 features and removed 59287 features. 

### Applying the filter

I am creating filtered versions of the counts matrix and features table using the filtering rule.

```{r}
counts_filtered <- counts[keep, ]
features_filtered <- features[keep, ]

dim(counts)
dim(counts_filtered)

dim(features)
dim(features_filtered)
```

After filtering, the counts matrix decreased from 89,118 to 29,831 feature rows, while all 181 sample columns remained. The features table also decreased to 29,831 rows while retaining its two columns. 

### Testing alternative thresholds

I am testing stricter filtering rules to determine how changing the count or sample requirement affects the number of retained features.

```{r}
keep_20_samples <- rowSums(counts >= 5) >= 20
keep_count_10 <- rowSums(counts >= 10) >= 10

c(
  original_rule = sum(keep),
  require_20_samples = sum(keep_20_samples),
  require_count_10 = sum(keep_count_10)
)
```

The original filtering rule: 29,831 features. 
At least five counts:26,607 features, 
At least ten counts retained 23,020 features. Both stricter rules removed more features. 

### Effect of filtering on feature types

I am comparing the numbers of canonical genes, LINE elements, and HERV elements before and after filtering.

```{r}
table(features$gene.type)
table(features_filtered$gene.type)
```
After filtering 25,050 canonical genes, 2,783 LINE elements, and 1,998 LTR/HERV elements remained. Canonical genes lost the largest number of features because this category originally contained the most However, LTR elements lost the largest proportion.

## Part 3 Checkpoint: What I Learned

Many features had very low expression across the entire dataset:

10% had zero total counts.
26% had fewer than 10 total counts.
50.25% had fewer than 50 total counts.

I used the filtering rule provided in the instructions. Keep a feature if it had at least five counts in at least ten samples.

The filter:

Kept 29,831 features
Removed 59,287 features
Retained 33% of the original matrix

The filtered counts matrix had 29,831 rows and 181 sample columns. The filtered features table also had 29,831 rows, confirming that the same rows were removed from both objects.

I also tested stricter rules:

Original rule: 29,831 features retained
At least five counts in 20 samples: 26,607 retained
At least ten counts in ten samples: 23,020 retained

Increasing the count requirement removed the most features. 

After filtering, the remaining feature types were:

25,050 canonical genes
2,783 LINE elements
1,998 LTR/HERV element

Canonical genes lost the greatest number of features because they were originally the largest category. However, LTR/HERV elements lost the greatest percentage.

## Part 4: Inferring sample sex from expression

Helena says I need to use XIST expression to look for two groups of samples. XIST is generally more highly expressed in samples with two X chromosomes and lower in samples with one X chromosome.

First we should check XIST and examine its values

```{r}
"XIST" %in% rownames(counts)

xist_counts <- counts["XIST", ]

summary(xist_counts)
sort(xist_counts)
```
Next need to create a bar plot

```{r}
barplot(
  sort(xist_counts),
  las = 2,
  cex.names = 0.3,
  main = "XIST expression (sorted)",
  ylab = "Raw counts",
  col = "steelblue"
)
```


I will also create a histogram

```{r}
hist(
  xist_counts,
  breaks = 30,
  main = "Distribution of XIST counts",
  xlab = "XIST raw counts",
  col = "steelblue"
)
```

XIST was present in the counts matrix. Its expression ranged from zero to 14171 counts. The sorted values showed a low expression group ending around 649 counts and a higher.

```{r}
counts[1:5, 100]
```

```{r}
counts[1:7, 100:101]
```

```{r}
counts[1:7, c(100,102)]
```

```{r}
xist_threshold <- 842

sex_assignment <- ifelse(
  xist_counts > xist_threshold,
  "Female",
  "Male"
)

table(sex_assignment)
```

```{r}
colors <- ifelse(
  sex_assignment == "Female",
  "red",
  "blue"
)

barplot(
  sort(xist_counts),
  las = 2,
  cex.names = 0.3,
  main = "XIST expression, colored by assigned sex",
  ylab = "Raw counts",
  col = colors[order(xist_counts)]
)

legend(
  "topleft",
  legend = c("Female", "Male"),
  fill = c("red", "blue")
)
```


```{r}
xist_threshold <- 78

sex_assignment <- ifelse(
  xist_counts > xist_threshold,
  "Female",
  "Male"
)

table(sex_assignment)
```

Mostly female. Need to ask Helena about the right threshold. ChatGPT said you need to do a midpoint between a big jump. But Perplexity AI and Claude AI stated that 78 was a fair threshold.

I next need to examine the Y-chromosome genes to evaluate whether the assignments were consistent.

```{r}
sex_table <- data.frame(
  sample = names(xist_counts),
  XIST_count = as.numeric(xist_counts),
  assigned_sex = unname(sex_assignment)
)

head(sex_table, 20)
```

```{r}
colors <- ifelse(
  sex_assignment == "Female",
  "red",
  "blue"
)

barplot(
  sort(xist_counts),
  las = 2,
  cex.names = 0.3,
  main = "XIST expression, colored by assigned sex",
  ylab = "Raw XIST counts",
  col = colors[order(xist_counts)]
)

legend(
  "topleft",
  legend = c("Female", "Male"),
  fill = c("red", "blue")
)
```

### Validating assignments with Y-chromosome genes

Samples assigned male should generally express Y-linked genes, while samples assigned female should have little or no Y-linked expression. I am checking five Y linked genes to determine whether they support the XIST based assignments.

```{r}
y_genes <- c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "UTY")

y_genes %in% rownames(counts)
```

```{r}
if ("RPS4Y1" %in% rownames(counts)) {
  
  rps4y1_counts <- counts["RPS4Y1", ]
  
  barplot(
    rps4y1_counts[order(xist_counts)],
    las = 2,
    cex.names = 0.3,
    main = "RPS4Y1 expression, sorted by XIST",
    ylab = "Raw counts",
    col = colors[order(xist_counts)]
  )
  
  legend(
    "topright",
    legend = c("Female", "Male"),
    fill = c("red", "blue")
  )
}
```
All five Y-chromosome genes were present in the matrix. RPS4Y1 was strongly expressed in nearly all samples assigned male. However, one female assigned sample had noticeable RPS4Y1 expression. We should test more genes.

```{r}
plot_y_gene <- function(gene) {
  
  gene_counts <- counts[gene, ]
  
  barplot(
    gene_counts[order(xist_counts)],
    las = 2,
    cex.names = 0.3,
    main = paste(gene, "expression, sorted by XIST"),
    ylab = "Raw counts",
    col = colors[order(xist_counts)]
  )
  
  legend(
    "topright",
    legend = c("Female", "Male"),
    fill = c("red", "blue")
  )
}
```

```{r}
plot_y_gene("EIF1AY")
plot_y_gene("DDX3Y")
plot_y_gene("KDM5D")
plot_y_gene("UTY")
```

```{r}
validation_table <- data.frame(
  sample = colnames(counts),
  XIST = as.numeric(xist_counts),
  assigned_sex = unname(sex_assignment),
  RPS4Y1 = as.numeric(counts["RPS4Y1", ]),
  EIF1AY = as.numeric(counts["EIF1AY", ]),
  DDX3Y = as.numeric(counts["DDX3Y", ]),
  KDM5D = as.numeric(counts["KDM5D", ]),
  UTY = as.numeric(counts["UTY", ])
)

validation_table$Y_total <- rowSums(validation_table[, 4:8])

female_validation <- validation_table[
  validation_table$assigned_sex == "Female",
]

head(
  female_validation[
    order(female_validation$Y_total, decreasing = TRUE),
  ],
  10
)
```

SCZ36 was the only disagreement. Although its XIST count of 636 caused it to be initially assigned female, it strongly expressed all five Y-linked genes with a combined total of 4,685 counts. The next highest Y gene total among female-assigned samples was only 56. 

```{r}
# Preserve the original XIST-only assignments
sex_assignment_xist <- sex_assignment

# Create the final validated assignments
sex_assignment_validated <- sex_assignment_xist
sex_assignment_validated["SCZ36_GSM7796477"] <- "Male"

table(sex_assignment_validated)
```

### Separating the matrix by sex assignment

I am separating the original counts matrix into female- and male-assigned samples using the XIST results validated with Y-chromosome genes.

```{r}
female_samples <- names(
  sex_assignment_validated[
    sex_assignment_validated == "Female"
  ]
)

male_samples <- names(
  sex_assignment_validated[
    sex_assignment_validated == "Male"
  ]
)

length(female_samples)
length(male_samples)
```

```{r}
counts_female <- counts[, female_samples]
counts_male <- counts[, male_samples]

dim(counts_female)
dim(counts_male)
```

The female matrix contained 89,118 features across 113 samples, while the male matrix contained the same 89,118 features across 68 samples. The two groups total 181 samples, confirming that every sample was included exactly once.

### Sex breakdown by condition

I am combining the validated sex assignments with the schizophrenia and control labels to determine the number of samples in each group.

```{r}
sex_condition <- data.frame(
  sample = colnames(counts),
  sex = unname(
    sex_assignment_validated[colnames(counts)]
  ),
  condition = ifelse(
    grepl("SCZ", colnames(counts)),
    "Schizophrenia",
    "Control"
  )
)

table(sex_condition$sex, sex_condition$condition)
```

```{r}
write.csv(
  sex_condition,
  "results/sex_assignments.csv",
  row.names = FALSE
)

saveRDS(
  counts_filtered,
  "results/dataset/counts_filtered.rds"
)

saveRDS(
  features_filtered,
  "results/dataset/features_filtered.rds"
)
```
```{r}
write.csv(
  sex_condition, file = "results/sex_assignments.csv",
  row.names = FALSE
)

saveRDS(
  counts_filtered,
  "results/dataset/counts_filtered.rds"
)

saveRDS(
  features_filtered,
  "results/dataset/features_filtered.rds"
)
```

## Part 4 Checkpoint: What I learned

The initial assignments were:

114 female
67 male

I then validated these assignments using five Y-chromosome genes:

All five genes were present in the matrix. Nearly all male-assigned samples strongly expressed the Y-linked genes, while nearly all female-assigned samples had little or no expression.

One sample, `SCZ36`, was initially assigned female because its XIST count was 636. However, it strongly expressed all five Y-linked genes, with a combined Y-gene count of 4,685. I  reclassified this sample as male.

The final validated assignments were:

113 female
68 male

I separated the counts matrix into:

A female matrix with 89,118 features and 113 samples
A male matrix with 89,118 features and 68 samples

The final breakdown by condition was:

Female Control - 67
Female SCZ - 46
Male Control - 30
Male SCZ - 38


## Part 6: Wrap Up

After speaking with Helena, I learned that R results can be displayed and saved in different ways. An R Notebook combines code, output, graphs, and written observations, while knitting turns the completed notebook into a shareable HTML file. We also discussed that sample sex should not be inferred from only one gene. I first used XIST to make initial assignments and then used five Y-chromosome markers to validate them. This showed the importance of using multiple methods to confirm a result and documenting the full process.


