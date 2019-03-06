# BSM research

## What's this?

* my computational research project on genomic imprinting in the brain
    * the body of this work is statistical analysis of RNA-seq data
* the code and documentation of the analysis a.k.a my **lab notebook**
* collection of scripts and small data sets for project-related computation
* slides for scientific presentations

## Introduction

This computational project is a large part of a biomedical research study at the Mount Sinai School of Medicine in Andy Chess' laboratory.  The study has been published as *Unperturbed expression bias of imprinted genes in schizophreniaa* ([Nat Commun. 2018 Jul 25;9(1):2914](https://www.nature.com/articles/s41467-018-04960-9)).

Along with the research article my lab notebook was published as a [website](https://attilagk.github.io/monoallelic-brain-notebook/), whose source is in my [monoallelic-brain-website](https://github.com/attilagk/monoallelic-brain-website) repository.  The present repository, however, contains not only the source of the articles of lab notebook but also scripts, or presentations (see next section) while it lacks the website specific files present in [monoallelic-brain-website](https://github.com/attilagk/monoallelic-brain-website).

## Structure and content

```
ms/             # manuscript of Nat Commun. 2018 Jul 25;9(1):2914
notebook/       # the main matter: my lab notebook
presentations/  # LaTeX/Beamer slides
src/            # scripts for project-wide use
```

As noted above the main matter of this repository is the lab notebook, a chronological sequence of articles in the `notebook/` directory.  Technically these articles are R Markdown documents, which means that they contain code chunks to be executed by interpreters such as `R`, or `bash`.  Longer, more complex code chunks have been moved to scripts located either in `notebook/` (these are associated to some article) or in `src/` (these are shared by multiple articles).  Finally, the most general-purpose code is part of [research-tools](https://github.com/attilagk/research-tools), a separate repository.

## The data

The [CommonMind Consortium](http://commonmind.org/) generates, stores and shares the RNA-seq and SNP-array data that have been analyzed during this work.
