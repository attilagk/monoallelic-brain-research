---
layout: post
---

## Introduction

Based on the aggregation of various p-values from [a previous analysis]({{ site.baseurl }}{% post_url 2016-10-03-permutation-test %}), a list of those genes was collated, which are significantly affected by any of the four chosen coefficients (Age, GenderMale, DxSCZ, Ancestry.1) representing various biological variables (see [signif-gene-effects-either.csv]).  These genes were then queried against Ensemble Genes 86 database, Homo sapiens genes (GRCh38.p7) using BioMart with the following HGNC symbols as input external reference IDs:

MAGEL2,
MIR668,
MEG8,
NAP1L5,
MEG3,
PEG3,
NDN,
PEG10,
KCNQ1OT1,
ZDBF2,
KCNK9,
INPP5F,
MEST,
PWRN1,
UBE3A

Note that in the above list MIR668 replaced AL132709.5 and MEG8 replaced RP11-909M7.3 because the latter IDs used as HGNC symbols gave no hit in Ensemble.

## Annotation

The following attributes were obtained:

Gene type,
Phenotype description,
Study External Reference,
Description,
Associated Gene Name,
Ensembl Family Description,
Chromosome Name,
Gene Start (bp)

The resulting table was saved in [signif-gene-effects-either-biomart.csv][signif-gene-effects-either-biomart.csv].  Those information were extended manually and saved in [signif-gene-effects-either-manual-annot.csv][signif-gene-effects-either-manual-annot.csv].

## Comparison with previous studies

### Copy number variation (CNV) in SCZ

Deletion or duplication of some 15 loci has been reported to be associated with SCZ ([Rees2014], [Sullivan2012]) but only 1 of these, in the 15q11q13 locus (Prader–Willi/Angelman region) contains a single gene whose parental bias is significantly affected by SCZ in our analysis: UBE3A.  The Prader–Willi/Angelman region overlaps the "SNRPN cluster" ([Peters2014]) of imprinted genes, whose dose variation, and in particular increased dose of the ubiquitin-protein ligase E3A (UBE3A), has been recognized as an important risk factor of psychosis and SCZ ([McNamara2013]).  The ubiquitination of synaptic proteins in particular GABAA has been found to be important for their turnover but UBE3A has also been linked to the glutamatergic system.

### Differential expression in SCZ

[Our comparison]({{ site.baseurl }}{% post_url 2016-10-20-differential-expression-scz %}) to Differentially Expressed Genes in Schizophrenia ([Fromer2016a]) shows that out of the 15 genes that our study found to be significantly associated to some biological predictor in terms of allelic bias, only one, PEG10, was found to be differentially expressed.  Moreover, none of the genes in the SNRPN cluster (i.e. the Prader-Willi/Angelman region) were among the $$\approx 700$$ differentially expressed genes.  This discrepancy might be resolved if, as suggested ([McNamara2013]), the dose of imprinted genes is more tightly regulated than that of non-imprinted genes, because that would mean that subtle perturbations in the expression of imprinted genes would be sufficient to cause phenotypic changes but insufficient to be detected in noisy RNA-seq measurements.  This further suggests that the read count ratio is a more sensitive indicator of altered expression of imprinted genes than the total read count.

### Animal models

We found MEST to be highly significantly associated to SCZ.  Disruption of MEST (a.k.a. PEG1) in mice showed involvement in embryonic growth as well as maternal behavior ([Lefebvre1998]).

[Rees2014]: http://www.ncbi.nlm.nih.gov/pubmed/24311552
[Sullivan2012]: http://www.ncbi.nlm.nih.gov/pubmed/22777127
[Peters2014]: https://www.ncbi.nlm.nih.gov/pubmed/24958438
[McNamara2013]: http://www.ncbi.nlm.nih.gov/pubmed/23697931
[Fromer2016a]: http://www.ncbi.nlm.nih.gov/pubmed/27668389
[Lefebvre1998]: http://www.ncbi.nlm.nih.gov/pubmed/9771709

[signif-gene-effects-either.csv]: {{ site.baseurl }}/assets/projects/monoallelic-brain/signif-gene-effects-either.csv
[signif-gene-effects-either-biomart.csv]: {{ site.baseurl }}/assets/projects/monoallelic-brain/signif-gene-effects-either-biomart.csv
[signif-gene-effects-either-manual-annot.csv]: {{ site.baseurl }}/assets/projects/monoallelic-brain/signif-gene-effects-either-manual-annot.csv
