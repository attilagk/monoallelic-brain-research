# More discussion on manuscript

## Overall

* What are the most interesting findings?

> Andy: (1) There are seem to be much fewer imprinted genes as some prev. studies suggested
> (Catherine Dulac et al). (2) Age dependence of imprinting.  (3) First(?) study using human
> samples.

* What needs further analysis?
    * regression on age of death: model selection?
    * Error rates for calling monoallelic expression necessary?

> We started discussing the necessity of error control for supporting the first conclusion above.
> I will write a more formal account on error control soon.

## Regression

$\mathrm{E}[ \mathrm{LOI\_R} ] = x_1 \beta_1 + ...$, where $x_1$ is the age of death and the null hypothesis is $\beta_1=0$, i.e. age has no impact on imprinting. $\mathrm{LOI\_R}$ seems to be a variable that in some way aggregates $\{S_g\}_g$ over 8 or 13 genes $g$.  But what is the definition of $\mathrm{LOI\_R}$?  What kind of aggregation is it (summation, pooling,...)?

> We looked at Ifat's R script for regression analysis and the definition of LOI_R.  When I get
> access to her files on the other server I'll look at them.

## Error rates

The manuscript provides no error rates for the classification of genes as mono or biallelically
expressed.

1. frequentist approach
    * $p$-values based on the null distribution of $S_g^{(i)}$
    * FDR control based on estimate of fraction of monoallelically expressed genes
2. Bayesian approach
    * probabilities of mono/biallelic expression: $\pi(m) + \pi(b) = 1$
    * prior prob. based on 
        a. estimate of fraction of monoallelically expressed genes
        b. distance from known imprinted genes $\pi(m | \mathrm{dist})$ (further extension: HMM)
    * posterior $\pi(m | x, y)$ given
        a. expression data $x$ (RNA-seq)
        b. genotype data $y$ (SNP-array)
        c. likelihood $f(x,y|m)$ for $m$ based on $x,y$
        d. prior $\pi(m)$

### Notes

In the frequentist approach we only need the likelihood function $f(x,y|b)$, whereas in the
Bayesian one we also need $f(x,y|m)$ (and the prior $\pi$, of course).

> Andy: permutation-derived null distribution of $S_g^{(i)}$ seems preferable instead of binomial assumption

The form of likelihood $f$ depends on the dependency structure of
the following variables:


### Error of genotype calling

A different kind of error rates **is** provided: error for calling genotypes (Figures
[error rate 1] and [error rate 2].

> Andy: discordant call is when RNA-seq suggests monoallelic expression and the Chip-array suggests
> heterozygosity

Figure [error rate 1].

* Does `AB` mean heterozygosity for a given SNP?
* `all calls`: within an individual?  For all individuals?
* How was `probability for AB` calculated?  What is `chip`?
* What are `err_AB_100` and `err_AB_7`?

> Andy: the acceptable minimum number of fragments covering a given SNP.

Figure [error rate 2].

* What are `data points`?
* What are numbers (1 to 100) next to plotted symbols?

## association of HLA genes to schizophrenia

Nonsignificant tendency for HLA-DQB1 was found.  Worth following up?

> Andy: not really worth it.

Imputation of HLA types uses two sources of info

1. SNPs (HIBAG)
2. RNA-seq (PHLAT)

[error rate 1]: http://katahdin.mssm.edu/ifat/web/cm/figures/error1.html
[error rate 2]: http://katahdin.mssm.edu/ifat/web/cm/figures/error2.html
