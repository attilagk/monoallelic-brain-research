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

### Already done

* Create table explaining the meaning of model names such as *wnlm.Q*
* Word carefully the causal interpretation: *allelic bias -> SCZ*

[ASEReadCounter]: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php
[Lappalainen 2013]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3918453/
[Hoffman 2016]: http://www.cell.com/cell-stem-cell/pdfExtended/S1934-5909(16)30401-5
[Castel 2015]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574606/
[Baran 2015]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4484390/
[Pirinen 2015]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4514921/
[MAMBA]: http://www.well.ox.ac.uk/~rivas/mamba/
[Knowles 2015]: http://biorxiv.org/content/early/2015/09/13/025874
