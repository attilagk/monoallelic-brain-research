---
layout: post
tags: [ literature ]
---

## Introduction

Based on the aggregation of various p-values from [a previous analysis]({{ site.baseurl }}{% post_url /monoallelic-brain/2016-10-03-permutation-test %}), a list of those genes was collated, which are significantly affected by any of the four chosen coefficients (Age, GenderMale, DxSCZ, Ancestry.1) representing various biological variables (see [signif-gene-effects-either.csv]).  These genes were then queried against Ensemble Genes 86 database, Homo sapiens genes (GRCh38.p7) using BioMart with the following HGNC symbols as input external reference IDs:

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

## General comparison with previous studies

### Copy number variation (CNV) in SCZ

Deletion or duplication of some 15 loci has been reported to be associated with SCZ ([Rees2014], [Sullivan2012]) but only 1 of these, in the 15q11q13 locus (Prader–Willi/Angelman region) contains a single gene whose parental bias is significantly affected by SCZ in our analysis: UBE3A.  The Prader–Willi/Angelman region overlaps the "SNRPN cluster" ([Peters2014]) of imprinted genes, whose dose variation, and in particular increased dose of the ubiquitin-protein ligase E3A (UBE3A), has been recognized as an important risk factor of psychosis and SCZ ([McNamara2013]).  The ubiquitination of synaptic proteins in particular GABAA has been found to be important for their turnover but UBE3A has also been linked to the glutamatergic system.

### Differential expression in SCZ

[Our comparison]({{ site.baseurl }}{% post_url /monoallelic-brain/2016-10-20-differential-expression-scz %}) to Differentially Expressed Genes in Schizophrenia ([Fromer2016a]) shows that out of the 15 genes that our study found to be significantly associated to some biological predictor in terms of allelic bias, only one, PEG10, was found to be differentially expressed.  Moreover, none of the genes in the SNRPN cluster (i.e. the Prader-Willi/Angelman region) were among the $$\approx 700$$ differentially expressed genes.  This discrepancy might be resolved if, as suggested ([McNamara2013]), the dose of imprinted genes is more tightly regulated than that of non-imprinted genes, because that would mean that subtle perturbations in the expression of imprinted genes would be sufficient to cause phenotypic changes but insufficient to be detected in noisy RNA-seq measurements.  This further suggests that the read count ratio is a more sensitive indicator of altered expression of imprinted genes than the total read count.

## Specific genes

### ZDBF2

Little is known about this zinc finger doman containing protein.  One study found  the long isoform of Zdbf2 (Liz) "potentially functions as both Zdbf2-coding RNA and cis-regulatory RNA" ([Duffie2014]).

### NAP1L5

Its functional role is not understood.

### PEG10

See the article *Deletion of Peg10, an imprinted gene acquired from a retrotransposon, causes early embryonic lethality* ([Ono2006]).

### MEST (PEG1)

We found MEST (a.k.a. PEG1) to be highly significantly associated to SCZ.  [Lefebvre1998] disrupted MEST in mouse embryonic stem cells and observed embryonic growth retardation as well as abnormal maternal behavior.  MEST has been implicated with the Silver-Russell syndrome ([OMIM 180860], see also [Peters2014]) based on maternal uniparental disomy for chromosome 7 (mUPD7) at the segment 7q31-qter in some patients with the syndrome.  The Silver-Russell syndrome is characterized by growth- but not mental retardation.

### KCNK9

The PubMed entries related to KCNK (Gene ID [51305](https://www.ncbi.nlm.nih.gov/gene/51305)) show that this gene (a.k.a. TASK3) has not only been linked to Birk-Barel mental retardation dysmorphism syndrome but also cancer.

### INPP5F

Little is known about the physiological role of this gene (inositol polyphosphate-5-phosphatase F) but its related PubMed entries (Gene ID [22876](https://www.ncbi.nlm.nih.gov/gene/22876)) show it has been linked to cancer.

### KCNQ1OT1

As the name shows (KCNQ1 opposite strand/antisense transcript 1), KCNQ1OT is on the opposite strand relative to KCNQ1.  The imprinting control region (ICR) of this cluster of genes (including KCNQ1, of course) is in an intron of KCNQ1 and the transcription of KCNQ1OT is a well-known mechanism of imprinting (the "ncRNA model of imprinting" in [Plasschaert2014]).  "Besides in Beckwith-Wiedemann syndrome, [...] the transcript also plays an important role in colorectal carcinogenesis" (from Entrez Gene, ID [10984](https://www.ncbi.nlm.nih.gov/gene/10984)).

### MEG3, RP11-909M7.3 (MEG8), AL132709.5 (MIR668)

These genes were determined in our study as part of the same cluster.

As for MEG3, Entrez Gene (ID [55384](https://www.ncbi.nlm.nih.gov/gene/55384)) notes that this lincRNA is a known tumor suppressor. RP11-909M7.3/MEG8 is a lincRNA with enigmatic physiological role.  See Entrez Gene (ID [79104](https://www.ncbi.nlm.nih.gov/gene/79104)).  AL132709.5 (MIR668), a microRNA, also has enigmatic physiological role.  See Entrez Gene (ID [768214](https://www.ncbi.nlm.nih.gov/gene/768214)).

### MAGEL2, NDN, PWRN1, UBE3A

These genes are all in one cluster in the Prader-Willi region.

[Peters2014] reviews the role of MAGEL2 and NDN in obesity.  Ndn mutant mice point to a role in neuronal development [Plasschaert2014] as its protein product "may suppress growth in postmitotic neurons" (Entrez Gene, ID [4692](https://www.ncbi.nlm.nih.gov/gene/4692).

Some roles of UBE3A in human neuropsychiological disorders has been mentioned above.  In addition to those, UBE3A is known to be linked to Angelman syndrome (reviewed by [Plasschaert2014]) and the Ube3a knockout mouse has been studied extensively and impairments have been found in synapse developement, long term potentiation, sleep and contextual learning ([Plasschaert2014] and [Peters2014]).

### PEG3

Peg3/Pw1 mutant mice show abnormal p53-mediated apoptosis [Broad2009].  More precisely, the normal sex-specific differences in apoptosis in some brain regions (mostly related to olfaction) were diminished by Peg3 inactivation.  This role of Peg3 in anatomical sexual dimorphism explains why it is required for normal maternal and sexual behavior.  See the general review [Plasschaert2014] and a more specific one by [Keverne2015], who therein also promotes the "coadaptation hypothesis" of imprinting.  That hypothesis has been challenged by [Haig2014].  Nonetheless, the putative role of PEG3 in sex-specific behavior is in line with our finding that PEG3's parental bias is significantly associated with gender.

[Rees2014]: http://www.ncbi.nlm.nih.gov/pubmed/24311552
[Sullivan2012]: http://www.ncbi.nlm.nih.gov/pubmed/22777127
[Peters2014]: https://www.ncbi.nlm.nih.gov/pubmed/24958438
[McNamara2013]: http://www.ncbi.nlm.nih.gov/pubmed/23697931
[Fromer2016a]: http://www.ncbi.nlm.nih.gov/pubmed/27668389
[Lefebvre1998]: http://www.ncbi.nlm.nih.gov/pubmed/9771709
[Keverne2015]: http://www.ncbi.nlm.nih.gov/pubmed/25404322
[OMIM 180860]: http://omim.org/entry/180860
[Broad2009]: http://www.ncbi.nlm.nih.gov/pubmed/19224563
[Plasschaert2014]: http://www.ncbi.nlm.nih.gov/pubmed/24757003
[Haig2014]: http://www.ncbi.nlm.nih.gov/pubmed/24129605
[Duffie2014]: http://www.ncbi.nlm.nih.gov/pubmed/24589776
[Ono2006]: http://www.ncbi.nlm.nih.gov/pubmed/16341224

[signif-gene-effects-either.csv]: {{ site.baseurl }}/assets/projects/monoallelic-brain/signif-gene-effects-either.csv
[signif-gene-effects-either-biomart.csv]: {{ site.baseurl }}/assets/projects/monoallelic-brain/signif-gene-effects-either-biomart.csv
[signif-gene-effects-either-manual-annot.csv]: {{ site.baseurl }}/assets/projects/monoallelic-brain/signif-gene-effects-either-manual-annot.csv
