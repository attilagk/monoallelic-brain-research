## Goal

It is crucial for publication to ensure that the genes selected for regression analysis were chosen on objective criteria.  This requires the reproduction of the set of genes that were called monoallelic in the [manuscript][ms] drafted by Ifat.  I'll refer to that in this document as *the previous manuscript*.  If perfect reproduction fails then minimum discordance is desired, and the regression analysis will **need to be repeated** accordingly with the new set of genes called monoallelic.

## Procedures from the previous manuscript

### Filters

The previous manuscript makes several statements on how the data were filtered before the gene were ranked and presented on [Fig 1][ifat fig 1]:

section *Methods: Allelic expression analyses*, p3
> there had to be $\ge 7$ reads, with a quality score $\ge 20$; a further gene-wide filter requires $20$ or more such reads for the gene to be assessable in a given individual

section *Methods: Allelic expression analyses*, p4
> Genes that were supported by fewer than 25 subjects were excluded.

section *Results: Pipeline*, p4
> there had to be $\ge 7$ reads; a further gene-wide filter requires 15 or more such reads for the gene to be assessable in a given individual

section *Results: Landscape of monoallelic expression[...]*, p6
> sufficient coverage defined as $\gt 20$ reads mapping to heterozygous SNPs per individual in at least ten individuals

In addition, for the regression analysis a stronger filtering was used:

section *Results: Relaxation of imprinting*, p8
> We examined genes where we had greater than 180 analyzable individuals and where 30% or more of those individuals displayed monoallelic expression defined by Sg >0.9. An analyzable individual was defined as one with at least 50 reads at one or more SNPs imputed to be heterozygous (Lhet cutoff of 0.95 as above).

#### My interpretation

Two kind of filters were used: I call them **read count-based** and **individual-based**.  These have the following properties:

1. the read count-based filter removes any such pair $(i,g)$ of individual $i$ and genes $g$ for which the total read count $n_{ig}<t_\mathrm{rc}$, where $t_\mathrm{rc}$ is what I call a **read count threshold**
   * Note that there is *another* read count-based filter, which removes pairs $(s,g)$ of SNP $s$ and gene $g$ in the same fashion as the previous one but with threshold 7.  However, I cannot investigate this filter until the SNP-wise read counts become easily accessible.
1. the individual-based filter removes any genes $g$ (across all individuals) if read count data involving $g$ are available on less than $t_\mathrm{ind}$ number of individuals

Note that in fact there are two read count-based filters: not only the one for individual SNPs and a

Given the above quotations from the previous manuscript the read count-based filter might have been applied with
$t_\mathrm{rc} = 20 \text{ or } 15 \text{ or } 21$, whereas the individual-based filter with $t_\mathrm{ind} = 25 \text{ or } 10$ in order to arrive at [Fig 1][ifat fig 1].

### Ranking and calling monoallelic expression

The **ranking** of genes is based on the following score. The score of gene $g$ is the fractions of individuals $i$ for whom $s_{ig} \gt 0.9$.  In other words the score of $g$ is based on the empirical cumulative distribution function ECDF $\hat{F}_g$ or, equivalently, on the survival function $1-\hat{F}_g$:
$$\text{the score of gene } g = 1 - \hat{F}_g(0.9).$$

Judged from [Fig 1][ifat fig 1] the **classification** of genes appears to have been defined as follows: $\text{the score of gene } g \ge 0.3 \Rightarrow g$ is called monoallelic.  Otherwise $g$ is called not monoallelic.  (Does "not monoallelic" imply biallelicity or do we wish to consider one or more intermediate classes?)  However, nothing is explicitly stated in the previous manuscript on the classification of genes independently of individuals.  For example (section *Results: Landscape of monoallelic expression[...]*, p6) we read
> 1b shows the 13 (of 32) genes that did not make it into the top Sg scoring group

but we are not told the criterion based on which the top scoring group was defined.

### Other quantities of interest

Although *not* used for ranking or classification, other fractions were also calculated and presented in [Fig 1][ifat fig 1].  First, not only $\hat{F}_g(s=0.9)$ was obtained but also $\hat{F}_g(s)$ for $s = 0.6,0.7,0.8$.  Second, the fraction of individuals were calculated that passes the test conceived by Andy:
$$ \text{given gene } g \text{ the fraction of indiv. passed} = \frac{\# \{i \,:\, s_{ig} \le 0.6 \text{ and } \mathrm{UCL}_{ig} \le 0.7 \}}{\# \{i\}}, $$
where the upper 95 % confidence limit is given by
$$\mathrm{UCL}_{ig} = s_{ig} + z_{0.975} \sqrt{\frac{s_{ig} (1 - s_{ig})}{n_{ig}}},$$
such that $z_{p}$ is the $p$ quantile of the standard normal distribution and $n_{ig}$ is the observed total read count.

## Results

### Genome-wide data import and preparation



Load functions:

```r
source("../../src/import-data.R")
source("../../src/utils.R")
```



Read count data have been imported (not shown).  We see that the minimum number of reads for any gene and individual is 7.

```r
min(unlist(N), na.rm = TRUE)
```

```
## [1] 7
```

The following code

1. applies the read count-based filter given a *sequence* of read count thresholds $t_\mathrm{rc}$
1. applies the individual-based filter given $t_\mathrm{ind}$
1. calculates the fractions of interest, which include the score $1-\hat{F}_g(0.9)$
1. ranks genes according to their score


```r
min.obs <- 25 # set t_ind
# implementation detail!: filter out genes with fewer observations than 'min.obs'
g.passed <- names(S)[sapply(S, function(y) sum(! is.na(y)) >= min.obs)]
# a sequence of thresholds t_rc for minimum read count to filter individual data points S_{ig} or N_{ig}
min.reads <- c(7, 15, 20)
min.r.names <- paste0("min.reads.", min.reads)
names(min.r.names) <- min.r.names
names(min.reads) <- min.r.names
# filter S_{ig} (then N_{ig}) first based on 'min.reads' and then again on 'min.obs'
Sf <- lapply(min.reads, function(m) filter.min.read(m, X = S[g.passed], N = N[g.passed], min.obs = min.obs))
Nf <- lapply(min.reads, function(m) filter.min.read(m, X = N[g.passed], N = N[g.passed], min.obs = min.obs))
# ECDFs for all filter levels and all genes g; individual ECDF components F_g are named according to gene g
ECDF <- lapply(Sf, sorted.ecdfs)
frac <- lapply(min.r.names,
               function(m)
                   do.fractions(ECDF[[m]], Sf[[m]], Nf[[m]],
                                frac = 10:6 / 10, ucl.fun = CI.p, max.ucl = 0.7, max.s = 0.6))
```
Note that the individual-based filter is also applied *before* the read count-based one.  However, this step is an implementation detail and has no conceptual significance.

#### Result on the individual-based filter

Besides $t_\mathrm{ind}=25$ the threshold $t_\mathrm{ind}=10$ was also tested (not shown) but the former resulted in clearly more consistency with Ifat's ranking (see below).

### Gene rankings

The plots show four gene rankings and the corresponding fractions of interest.  The first three correspond to the sequence of three $t_\mathrm{rc}$ settings of the read count-based filter so these rankings will be named $R_{7}, R_{15}, R_{20}$.  The fourth plot follows the ranking seen on Ifat's [Fig 1][ifat fig 1], and this ranking will be referred to as $R_\mathrm{Ifat}$.  Since that figure shows only the top 51 genes, the same is done here for also the first three plots.  Several genes in the top 51 according to $R_{7}, R_{15}$ or $R_{20}$ are missing from Ifat's top 51 and therefore the imprinting status for those genes will need to be obtained from elsewhere.  Given this lack of information the first three graphs are colored with orange.



![plot of chunk compare-to-ifats-fig](figure/compare-to-ifats-fig-1.png)

#### Gene rankings with Ifat's ranking as reference


```
##               imprinting.status R.ifat  R.7 R.15 R.20
## MAGEL2                    known      1   19    1    1
## TMEM261P1       candidate, <1MB      2    3    2   NA
## SNHG14                    known      3    1    3    2
## AL132709.5      candidate, <1MB      4    6    4    4
## RP11-909M7.3    candidate, <1MB      5   13    5    3
## NAP1L5                    known      6   10    7    6
## ZIM2                      known      7    2    6    5
## MEG3                      known      8    4    8    7
## PEG3                      known      9    5    9    8
## PWAR6           candidate, <1MB     10    7   10    9
## FAM50B                    known     11   16   11   13
## NDN                       known     12    8   12   10
## SNURF                     known     13    9   13   11
## PEG10                     known     14   11   14   14
## SNRPN                     known     15   12   15   12
## KCNQ1OT1                  known     16   21   16   15
## ZDBF2                     known     17   15   17   16
## GRB10                     known     18   24   18   17
## SNORD116-20     candidate, <1MB     19   22   19   18
## KCNK9                     known     20   28   20   20
## INPP5F                    known     21   25   21   19
## HLA-DRB5              candidate     22   30   22   21
## RP13-487P22.1   candidate, <1MB     23   39   23   22
## GSTM1                 candidate     24   32   24   24
## MEST                      known     25   38   26   23
## hsa-mir-335     candidate, <1MB     26   56   31   NA
## IL1RL1                candidate     27   49   28   NA
## ZNF331                    known     28   41   29   25
## DIRAS3                    known     29   44   32   NA
## PWRN1           candidate, <1MB     30   43   33   29
## HLA-DQB1              candidate     31   50   35   31
## PAX8-AS1              candidate     32   52   37   33
## HNRNPU                candidate     33   53   38   32
## HLA-DQA1              candidate     34   65   44   47
## RP11-54F2.1           candidate     35   83   48   40
## SYT7                  candidate     36   81   46   34
## NME1-NME2             candidate     37   64   40   35
## RAD23A                candidate     38   95   57   NA
## NLRP2                     known     39  122   77   66
## IGF2                      known     40   84   83   82
## UBE3A                     known     41  132  112   79
## NTM                       known     42  318  177  150
## DGCR6                     known     43 1683 1310   NA
## OSBPL5                    known     44 2479  946  659
## NAA60                     known     45 1192 2251 2345
## DGCR6L                    known     46 3483 2392 1890
## BEGAIN                    known     47 5004 2448 2397
## AIM1                      known     48 2028 4490   NA
## DLGAP2                    known     49 3041 2382 2104
## GNAS                      known     50 5529 3754 2997
## ZFAT                      known     51 4714 5016 4399
```

Clearly, when the read count-based filter is used at $t_\mathrm{rc}=15$ there is a close match with Ifat's ranking.  In fact, the two rankings agree for the top 24 genes (except that ranks 6 and 7 are swapped):

```r
all.equal(genes.ifat.ranks[ , "R.15"][c(1:5, 7:6, 8:24)], 1:24)
```

```
## [1] TRUE
```

From this result we may conclude that Ifat's filter settings were $t_\mathrm{rc}=15$ and $t_\mathrm{ind}=25$ and the discrepancies we see might partly or entirely be due to rounding errors introduced at the export of $S$ values to csv files and/or different implementation of the calculation of gene scores.

### Discrepant genes

From now on we only scrutinize the ranking obtained with $t_\mathrm{rc}=15$ i.e.$R_{15}$ and compare it to Ifat's ranking i.e. $R_\mathrm{Ifat}$.  On what genes do these two rankings disagree?  In particular, with what genes, if any, should I extend my latest (already extended) regression analysis?

Notice that the lowest scoring gene in the "candidate, < 1MB" category is *PWRN1* (green on [Fig 1][ifat fig 1]).  Its score is 0.5 and it ranks at 33 according to $R_{15}$ and at `genes.ifat.ranks["PWRN1", "R.ifat"]`.  We may define the classification threshold such that we only call genes monoallelic if their $R_{15}$-ranking is above *PWRN1*.

```r
d.genes <- list() # discrepant genes
(d.genes$A <-
    names(ECDF$min.reads.15)[ setdiff(seq_len(genes.ifat.ranks["PWRN1", "R.15"]),
                                      genes.ifat.ranks$R.15[seq_len(genes.ifat.ranks["PWRN1", "R.ifat"])]) ])
```

```
## [1] "RPL23AP7" "MRPL28"   "MRPS34"
```

However, adjusting the classification threshold to *PWRN1* looks quite arbitrary.  If we instead call genes with score $\gt 0.5$ monoallelic then one more gene, *GFRA2*, must also called:

```r
(d.genes$B <-
    setdiff(names(ECDF$min.reads.15)[sapply(ECDF$min.reads.15, function(f) 1 - f(0.9) >= 0.5) ],
            names(ECDF$min.reads.15)[ seq_len(match("PWRN1", names(ECDF$min.reads.15))) ]))
```

```
## [1] "GFRA2"
```

If we further lower the classification threshold to $0.3$ then according to $R_{15}$ these genes must also be called monoallelic.

```r
(d.genes$C <-
    setdiff(names(ECDF$min.reads.15)[sapply(ECDF$min.reads.15, function(f) 1 - f(0.9) >= 0.3) ],
            names(ECDF$min.reads.15)[sapply(ECDF$min.reads.15, function(f) 1 - f(0.9) >= 0.5) ]))
```

```
##  [1] "HLA-DQB1"      "HLA-DRB1"      "PAX8-AS1"      "HNRNPU"       
##  [5] "ZDHHC11"       "NME1-NME2"     "AHNAK"         "PPP1R37"      
##  [9] "PRKAR1B"       "HLA-DQA1"      "ANKRD18A"      "SYT7"         
## [13] "TMIE"          "RP11-54F2.1"   "MACF1"         "GPR75-ASB3"   
## [17] "THOC3"         "AC008443.1"    "TRIM52"        "CTXN1"        
## [21] "RP11-215G15.5" "UBAC1"         "RAD23A"
```

Some of these genes like *HLA-DQB1* are present in the top 51 according to the ranking $R_\mathrm{ifat}$ and hence their imprinting status can be read off from that table. This shows that all of those belong to the "candidate" category.  Looking at the status of all discrepant genes established above, with the most liberal classification threshold 0.3, we see that all of them are "candidate" (red) and therefore need not be included in further regression analysis.

```r
levels(gene.summary$imprinted) <- rev(c("known", "candidate, <1MB", "candidate")) # follow new nomenclature
d.genes.imprint <- as.character(gene.summary[unlist(d.genes), "imprinted"])
names(d.genes.imprint) <- unlist(d.genes)
d.genes.imprint
```

```
##      RPL23AP7        MRPL28        MRPS34         GFRA2      HLA-DQB1 
##   "candidate"   "candidate"   "candidate"   "candidate"   "candidate" 
##      HLA-DRB1      PAX8-AS1        HNRNPU       ZDHHC11     NME1-NME2 
##   "candidate"   "candidate"   "candidate"   "candidate"   "candidate" 
##         AHNAK       PPP1R37       PRKAR1B      HLA-DQA1      ANKRD18A 
##   "candidate"   "candidate"   "candidate"   "candidate"   "candidate" 
##          SYT7          TMIE   RP11-54F2.1         MACF1    GPR75-ASB3 
##   "candidate"   "candidate"   "candidate"   "candidate"   "candidate" 
##         THOC3    AC008443.1        TRIM52         CTXN1 RP11-215G15.5 
##   "candidate"   "candidate"   "candidate"   "candidate"   "candidate" 
##         UBAC1        RAD23A 
##   "candidate"   "candidate"
```

## Conclusion

* Ifat's ranking is the most consistent with the present $R_{15}$ ranking obtained with filter settings were $t_\mathrm{rc}=15$ and $t_\mathrm{ind}=25$
* the consistency is not perfect, which might be due to implementation details
* therefore the two ranking leads to discrepant sets of monoallelically called genes
* the level of discrepancy depends on what rule is chosen for classification
* even with the most liberal rule (the one yielding the largest called monoallelically set), all called genes fall in the "candidate" category
* thus the **regression** analysis **need not be further extended**


[ifat fig 1]: https://docs.google.com/presentation/d/1YvpA1AJ-zzir1Iw0F25tO9x8gkSAzqaO4fjB7K3zBhE/edit#slide=id.p4
[ms]: https://docs.google.com/document/d/1cWd4UH98SJR5lihDihC0ZO-C_A1-8MQ5COcixxCLzHE/edit
