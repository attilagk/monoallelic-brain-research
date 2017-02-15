---
layout: post
tags: [ andy, manuscript, gabrielhoffman, discussion ]
---

1. Extend regression analysis (both fixed and mixed effects models) to many more (perhaps all?) assessable genes.
1. Consolidate and employ both fixed and mixed effects models
    * $$\hat{\beta}$$ from mixed models to $$\hat{\beta}$$ from the corresponding fixed models in order to resolve present discrepancies.
    * use `varPartConfInf` for confidence intervals
1. Include HPCC data set; Gabriel didn't think it was as crucial as the improvements in statistical methodology he suggested.
1. Adapt filtering/correction procedures for cleaner read counts
    * consult with papers from Lappalainen et al: [Lappalainen 2013 Nature][Lappalainen 2013], [Castel 2015][Castel 2015] (tools and best practices), [Baran 2015][Baran 2015] (the landscape of genomic imprinting in human tissues) as well as Gabriel's methodological description ([Hoffman 2016][Hoffman 2016]) based on the above mentioned studies
    * reference allele bias
    * imputed dosage between 0.99 and 1.01
1. Extend predictors (covariates) with additional ones present in CMC
1. Fit simple regression models akin to Figure 5 "Fitting four families of generalized linear models..." in order to demonstrate qualitatively the results from multiple regression models.  In other words, extend Fig 5 to more genes and predictors.
1. Consider other tools and statistical approaches
    * the GATK walker [ASEReadCounter][ASEReadCounter]
    * the Bayesian model comparison framework by [Pirinen 2015] and the [MAMBA] software implementing it
    * [Knowles 2015][Knowles 2015]
1. Somehow use available knowledge on DNA methylation and age.  A name (Steve Hoggarth?) was mentioned.

## Already done

* Create table explaining the meaning of model names such as *wnlm.Q*
* Word carefully the causal interpretation: *allelic bias -> SCZ*

## Possible improvements

Here I assume we want to keep our earlier two-step approach:

1. obtain higher and lower read counts $$H_{ig}, L_{ig}$$
    * correct for various errors
    * aggregate over all het-SNP sites within gene $$g$$
1. modeling $$H_{ig}, L_{ig}$$ in order to call imprinted genes and study dependence on predictors 

Note that a few one-step approaches have been published, in which a single statistical model explicitly incorporates all sources of technical variation (errors), biological variation, and performs aggregation.  These one-step approaches are more powerful and robust but are not readily applicable to our specific dataset and aims.

### Obtaining higher and lower read counts $$H_{ig}, L_{ig}$$

* reference allele bias (mapping)
* genotyping/impputation error
    * previously used by us: a qualitative post-hoc test (ref/non-ref test) that is subjective therefore not implemented computationally
    * the fact that such post-hoc correction was necessary suggests that genotyping/imputation error was not sufficiently taken into account
* phasing error (aggregation)
    * the assignment of each SNP variant $$v_s \in \{ \mathrm{ref}_s, \mathrm{nonref}_s \}$$ at site $$s$$ to the higher (lower) expressed allele $$\{a_\mathrm{H}, a_\mathrm{L}\}$$ based on the higher and lower read counts at $$s$$: $$H_s, L_s$$.  Our current simple rule is that if $$\#\{\mathrm{ref}_s\} = H_s$$ then $$\mathrm{ref}_s = a_\mathrm{H}$$ otherwise $$\mathrm{ref}_s = a_\mathrm{L}$$

For some of these points there exist recommended tools and best practices (see [Castel 2015] and references therein).  For the phasing error various approximations can be envisioned based on probabilistic sampling from the space of haplotypes given some estimate of the expected allele ratio for transcripts.

### Modeling $$H_{ig}, L_{ig}$$

Previously we considered the ratio $$S_{ig} = H_{ig} / T_{ig}$$ where $$T_{ig} = H_{ig} + L_{ig})$$.  It was challenging to find a generalized linear regression model (GLM) that fits data $$\{S_{ig}, T_{ig} \, : \, i = 1,...,I \}$$ for all genes $$g$$, where $I$ is the number of individuals (samples).  No single link function (e.g. linear, sigmoidal or scaled sigmoidal) appeared good enough and no error distribution could capture the **overdispersion**: the binomial distribution of the logi.S and logi2.S was of course better than the homoscedastic normal model, but not good enough in general.  Eventually, a quasi-log transformation was applied to $$S_{ig}$$ and normal linear model was fitted (wnlm.Q) based on model checking diagnostics.

Two alternatives are sketched below, both of which models $$L_{ig}$$ given $$H_{ig}$$ instead of $$S_{ig}$$ given $$T_{ig}$$.  Both alternatives are inspired by published and well-tested statistical approaches and therefore

* both approaches might afford more power and robustness than our previous approach using the wnlm.Q and logi.S models
* they might be general enough to fit data for **all genes** (assessable) and not only the earlier selected 30 genes
* therefore, they may allow both genome-wide ranking of genes and genome-wide investigation of how allelic bias depends on biological predictors

#### Approach 1: negative binomial regression

The defining features are

1. log link function: $$\mu_{ig} = e^{x_{i}^\top \beta_g}$$
1. negative binomial distribution of $$L_{ig}$$ given its mean $$\mu_{ig}$$ as well as parameter $$\nu_{ig}$$, which controls dispersion.

$$
\begin{equation}
f(L_{ig}; \mu_{ig}, \nu_{ig}) = \frac{\Gamma(L_{ig} + \nu_{ig})}{\Gamma(\nu_{ig}) L_{ig}!}
\frac{\mu_{ig}^{L_{ig}} \nu_{ig}^{\nu_{ig}}}{(\mu_{ig} + \nu_{ig})^{L_{ig} + \nu_{ig}}}
\end{equation}
$$

Unfortunately this model cannot be readily employed to our data because the lower read counts $$L_{ig}$$ also depend on $$H_{ig}$$ (or equivalently on $$T_{ig}$$).  Perhaps $$L_{ig}$$ could be normalized using $$H_{ig}$$ or $$T_{ig}$$ separately for each gene $$g$$.

Note that parameter $$L_{ig}$$ may be interpreted as the number of failures whereas $$\nu_{ig}$$ the number of successes.  But $$H_{ig}$$ might also be given the same interpretation (number of successes).  What is the relationship between $$H_{ig}$$ and $$\nu_{ig}$$?  It's something to unfold.

#### Approach 2: transformation and linear model

1. linear link function: $$\mu_{ig} = x_{i}^\top \beta_g$$
1. normal distribution of the log odds $$\log \frac{L_{ig}}{H_{ig}}$$ given its mean $$\mu_{ig}$$ and variance $$\sigma^2_g / w_{ig}$$, where $$w_{ig}$$ is precision weight.

Challenge: estimation of weights $$w_{ig}$$; (how) could be done using voom/limma?



[ASEReadCounter]: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php
[Lappalainen 2013]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3918453/
[Hoffman 2016]: http://www.cell.com/cell-stem-cell/pdfExtended/S1934-5909(16)30401-5
[Castel 2015]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574606/
[Baran 2015]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484390/
[Pirinen 2015]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514921/
[MAMBA]: http://www.well.ox.ac.uk/~rivas/mamba/
[Knowles 2015]: http://biorxiv.org/content/early/2015/09/13/025874
