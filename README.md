---
title: "README"
author: "Jose V. Die"
date: "10/25/2024"
output: 
  html_document:
    keep_md: yes  
---



Load libraries

```r
library(Biostrings)
library(stringr)
library(dplyr)
library(tibble)
```

Load FASTA file as DNAStringSet object

```r
p = readDNAStringSet("protein.fasta")
tmp = names(p)
```

Look at the FASTA file 

```
## DNAStringSet object of length 35836:
##         width seq                                           names               
##     [1]   353 MARSGDTVGSKGVKRSKTRAY...RWGDNRDRSSNDTTHYTNGRT GWHPBRAD000001	mR...
##     [2]   353 MARSGDTVGSKGVKRSKTRAY...RWGDNRDRSSNDTTHYTNGRT GWHPBRAD000002	mR...
##     [3]   163 MDSSGGNCNVKKYKVGDGCGG...KNVRRKNNRSRHVKVDSGYWS GWHPBRAD000003	mR...
##     [4]   224 MRNRSSRVGKGSRNKMAATAS...SAGDMDSKMVDGSDSDSSSSA GWHPBRAD000004	mR...
##     [5]    98 MGVSAAAAAVDSGSAVTSDSH...ASKDDVKKKTYKRWAKKKSKT GWHPBRAD000005	mR...
##     ...   ... ...
## [35832]   574 MATGASVKVKRHGTAVTSAAC...AADHDHSNSNGDGVAHRKGRC GWHPBRAD035832	mR...
## [35833]    69 MVVDGSKDDTTATSGCCSNAT...TKADYAKDGNVSSSSSMSSTT GWHPBRAD035833	mR...
## [35834]   338 MGDRWVNDRAMVSVGTACVSS...TNDKVVHTSTGGAAKKVGASK GWHPBRAD035834	mR...
## [35835]   279 MTSSASSSAVVVSTGKDVAAS...HKKKTWRTTKSHDNTAVMCSR GWHPBRAD035835	mR...
## [35836]   572 MTDKNGSRATTDMWTSKWNDR...HKKKTWRTTKSHDNTAVMCSR GWHPBRAD035836	mR...
```



Clean protein names

```r
tmp = names(p)
tmp = str_split(str_squish(tmp), "Gene=")
```


```r
tmp = sapply(1:length(tmp), function(i) str_split_1(tmp[[i]][2], " ")[1])
names(p) = tmp
p
```

```
## DNAStringSet object of length 35836:
##         width seq                                           names               
##     [1]   353 MARSGDTVGSKGVKRSKTRAY...RWGDNRDRSSNDTTHYTNGRT GWHGBRAD000001
##     [2]   353 MARSGDTVGSKGVKRSKTRAY...RWGDNRDRSSNDTTHYTNGRT GWHGBRAD000001
##     [3]   163 MDSSGGNCNVKKYKVGDGCGG...KNVRRKNNRSRHVKVDSGYWS GWHGBRAD000002
##     [4]   224 MRNRSSRVGKGSRNKMAATAS...SAGDMDSKMVDGSDSDSSSSA GWHGBRAD000003
##     [5]    98 MGVSAAAAAVDSGSAVTSDSH...ASKDDVKKKTYKRWAKKKSKT GWHGBRAD000004
##     ...   ... ...
## [35832]   574 MATGASVKVKRHGTAVTSAAC...AADHDHSNSNGDGVAHRKGRC GWHGBRAD028623
## [35833]    69 MVVDGSKDDTTATSGCCSNAT...TKADYAKDGNVSSSSSMSSTT GWHGBRAD028624
## [35834]   338 MGDRWVNDRAMVSVGTACVSS...TNDKVVHTSTGGAAKKVGASK GWHGBRAD028625
## [35835]   279 MTSSASSSAVVVSTGKDVAAS...HKKKTWRTTKSHDNTAVMCSR GWHGBRAD028626
## [35836]   572 MTDKNGSRATTDMWTSKWNDR...HKKKTWRTTKSHDNTAVMCSR GWHGBRAD028626
```



```r
# Build a dataframe with n. isoforms and longest protein size per gene
df = tibble::tibble(gene= names(p), aa= width(p))

df %>% 
  arrange(gene, desc(aa)) %>% 
  count(gene, sort = T) %>% 
  filter(n > 1)
```

```
## # A tibble: 4,623 × 2
##    gene               n
##    <chr>          <int>
##  1 GWHGBRAD003538    15
##  2 GWHGBRAD006368    15
##  3 GWHGBRAD011440    14
##  4 GWHGBRAD005819    13
##  5 GWHGBRAD006276    12
##  6 GWHGBRAD012733    12
##  7 GWHGBRAD012738    11
##  8 GWHGBRAD017471    11
##  9 GWHGBRAD001874     9
## 10 GWHGBRAD005623     9
## # ℹ 4,613 more rows
```

```r
nisof <- df %>% 
  arrange(gene, desc(aa)) %>% 
  count(gene) %>% 
  pull(n)

df <- df %>%  
  arrange(gene, desc(aa)) %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  mutate(n_isof = nisof)

df %>% arrange(desc(n_isof))
```

```
## # A tibble: 28,626 × 3
##    gene              aa n_isof
##    <chr>          <int>  <int>
##  1 GWHGBRAD003538   443     15
##  2 GWHGBRAD006368   289     15
##  3 GWHGBRAD011440   682     14
##  4 GWHGBRAD005819   536     13
##  5 GWHGBRAD006276   288     12
##  6 GWHGBRAD012733   503     12
##  7 GWHGBRAD012738   658     11
##  8 GWHGBRAD017471   779     11
##  9 GWHGBRAD001874   638      9
## 10 GWHGBRAD005623   555      9
## # ℹ 28,616 more rows
```


```r
# Split sequences by their names
grouped_seqs <- split(p, names(p))
```


```r
# For each group, find the longest sequence
longest_seqs <- lapply(grouped_seqs, function(group) {
  group[which.max(width(group))]
})
```


```r
# Unlist the result to get back a DNAStringSet
longest_dna_set <- do.call(c, longest_seqs)
```


```r
# Assign correct names back to the sequences
names(longest_dna_set) <- names(longest_seqs)
longest_dna_set %>% head(3)
```

```
## $GWHGBRAD000001
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]   353 MARSGDTVGSKGVKRSKTRAYSG...CTRWGDNRDRSSNDTTHYTNGRT GWHGBRAD000001
## 
## $GWHGBRAD000002
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]   163 MDSSGGNCNVKKYKVGDGCGGST...VKKNVRRKNNRSRHVKVDSGYWS GWHGBRAD000002
## 
## $GWHGBRAD000003
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]   224 MRNRSSRVGKGSRNKMAATASTC...HTSAGDMDSKMVDGSDSDSSSSA GWHGBRAD000003
```



```r
# Convert the list back into a DNAStringSet
LargeDNAStringSet <- DNAStringSet(sapply(longest_dna_set, `[[`, 1))
```



```r
LargeDNAStringSet %>% head(3)
```

```
## DNAStringSet object of length 3:
##     width seq                                               names               
## [1]   353 MARSGDTVGSKGVKRSKTRAYSG...CTRWGDNRDRSSNDTTHYTNGRT GWHGBRAD000001
## [2]   163 MDSSGGNCNVKKYKVGDGCGGST...VKKNVRRKNNRSRHVKVDSGYWS GWHGBRAD000002
## [3]   224 MRNRSSRVGKGSRNKMAATASTC...HTSAGDMDSKMVDGSDSDSSSSA GWHGBRAD000003
```


```r
# Save a FASTA file 
writeXStringSet(x = LargeDNAStringSet, filepath = "protein_cnames.fa")
```



```r
cnames = readDNAStringSet("protein_cnames.fa")
cnames["GWHGBRAD003538"]
```

```
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]   443 MTSTAVWHMKGSKHRSRAYCDWC...KRATTDANDSVDMGMVGMGRRDC GWHGBRAD003538
```


```r
cnames["GWHGBRAD028626"]
```

```
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]   572 MTDKNGSRATTDMWTSKWNDRKT...YYHKKKTWRTTKSHDNTAVMCSR GWHGBRAD028626
```

