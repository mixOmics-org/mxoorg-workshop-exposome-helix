# Workshop: multivariate analysis of the exposome

### Author: Prof Kim-Anh Lê Cao

A hands-on practical applying multivariate analysis to exposome and multi-omics data
with the [mixOmics](https://mixomics.org) R package, built for the ISGlobal exposome
workshop. The session runs for about 1.5 hours and is built around a single practical
that participants work through on their own laptops, using data from the HELIX study.

| Audience | Prerequisites | Duration |
| --- | --- | --- |
| Biologists and epidemiologists | A working knowledge of R: you should be comfortable running code and reading plots. No programming is required | ~ 1.5 hours |

### Material

- [**Practical link**](https://guides.mixomics.org/mxoorg-workshop-exposome-helix/practical/exposome_analysis.html)

- [**Data**](Data/) — one `.RData` file, and the script that builds it

### What the session covers

The practical opens with a short quick start, running PCA to establish the workflow
of every method in the package: run the method, plot the samples, plot the variables.
It then works through three case studies, all on the same two data sets so that only
the method changes:

1. **PCA** on the HELIX proteome and on the organochlorine exposures — unsupervised
   exploration of each in turn, colouring the children by sex, recruitment cohort and
   BMI category to help interpret the result.
2. **PLS1** — regression of a single exposure, summed serum PCBs, on the proteome,
   then sparse PLS1 to select the proteins that carry the association.
3. **block sPLS** — the same exposure family against four molecular blocks at once
   (proteome, serum metabolome, urine metabolome, transcriptome), the multiblock
   regression counterpart of DIABLO.

A **DIABLO** analysis, the classification counterpart of the third case study, can be
found in the companion practical at
[mxoorg-workshop-intro-mixomics](https://github.com/mixOmics-org/mxoorg-workshop-intro-mixomics).

### Before the workshop

**1. Install the software.** Install R, then RStudio. Use recent versions of both:

- [R](https://cran.r-project.org/) (R 4.0 or later)
- [RStudio](https://posit.co/download/rstudio-desktop/#download)

Then install mixOmics from Bioconductor, and check that it loads by typing the
following in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("mixOmics")

# No error messages means you are set. Warnings are fine.
library(mixOmics)
```

**Apple mac users:** if the imported `rgl` package will not install, install
[XQuartz](https://www.xquartz.org) first, then try again.

**2. Download this repository, and load the data once.** Clone it, or use
**Code > Download ZIP** on GitHub. The data file is about 7 MB, so please do this
**before** the session rather than over the venue's wifi at the start. Then check it
loads:

```r
# with your working directory set to the practical/ folder
load(file.path("..", "Data", "helix.RData"))
names(helix)
```

You should see `omics`, `exposure`, `outcomes`, `covariates`, `annotation` and
`provenance`.

---

## The data

### Where it comes from

The data come from **HELIX** (Human Early-Life Exposome), released by ISGlobal for its
**Exposome Data Challenge 2021**. HELIX is a collaborative project across six
established longitudinal population-based birth cohorts in six European countries:
France, Greece, Lithuania, Norway, Spain and the United Kingdom. Mother–child pairs
were followed from pregnancy, and the children were re-examined at roughly age 6 to 11.

What makes the dataset unusual is its design. Each child is characterised broadly —
across chemical, outdoor, indoor and lifestyle exposures — with four molecular layers
measured on the same children. That is the exposome idea: rather than one exposure and
one disease, the whole measured environment at once.

The extract distributed here is **not** the underlying cohort data, which is available
only on request subject to ethical and legislative review. See the
[HELIX data inventory](https://www.projecthelix.eu/index.php/es/data-inventory).
Investigators planning analyses similar to HELIX's own are encouraged to contact the
consortium at `helixdata@isglobal.org`.

### What is in the extract

`Data/helix.RData` contains a single list, `helix`, with 896 children — the
intersection of all five blocks, where the transcriptome is the binding constraint.
All 896 are in the **same order in every element**.

```
helix
├── omics           # four matrices, samples in rows
│   ├── proteome        896 x 36      circulating proteins, largely inflammatory markers
│   ├── serum           896 x 177     serum metabolites
│   ├── urine           896 x 44      urine metabolites
│   └── transcriptome   896 x 1000    transcripts (reduced, see below)
├── exposure        # the nine postnatal families, 81 variables, all numeric
│   ├── metals                896 x 9
│   ├── organochlorines       896 x 9
│   ├── phthalates            896 x 11
│   ├── phenols               896 x 7
│   ├── pfas                  896 x 5
│   ├── op_pesticides         896 x 4
│   ├── air_pollution         896 x 12
│   ├── built_environment     896 x 15
│   └── meteorological        896 x 9
├── outcomes        # data frame, all factors -- categorical
│   ├── sex          female 421, male 475
│   ├── cohort       6 levels: 152 / 83 / 165 / 147 / 199 / 150
│   ├── bmi_cat      thinness 9, normal 630, overweight 165, obese 92
│   ├── asthma       no 804, yes 92
│   └── mat_edu      primary 135, secondary 308, university 453
├── covariates      # data frame, numeric
│   └── child_age, zbmi_who, birthweight, raven_score, behaviour_total
├── annotation      # feature tables, each keyed by its block's column names
│   ├── proteome        Prot_ID, UniProtKB, Gene_Symbol, Gene_Name
│   ├── serum           Class, Class_2
│   ├── urine           var, CHEBI, KEGG
│   ├── transcriptome   GeneSymbol_Affy, EntrezeGeneID_Affy, CallRate
│   └── codebook        one row per exposure variable, with a description
└── provenance      # source, date, required citation, every processing step
```

### What was and was not done to it

**No block was transformed here.** All four arrived already transformed, and in three
different ways: the proteome and serum metabolome are log-transformed, the urine
metabolome is log2 with a quantification floor at 0.05 µM, and the transcriptome is
log2 and **per-gene centred**. Those states were established from the data rather than
read off the upstream documentation, which is silent on urine.

**No filtering removed anything except in the transcriptome.** The release ships
complete, imputed tables, and the three small blocks are curated panels of fewer than
200 features.

**The transcriptome was reduced for size, not for science.** It ships at 28,738
transcripts, was reduced to 5,000 in the scoping work (call rate ≥ 80%, then variance),
and is reduced here to the **top 1,000 of those by variance**. The 1,000 retain 55.6%
of the 5,000-probe variance and every one of them has a gene symbol, so loadings stay
interpretable. The reason is download size: the full 5,000 makes this file about 34 MB
against roughly 7 MB for 1,000, which matters when a room of people clones the
repository at once. This is a teaching-material decision and should not be mistaken
for a scientific one.

**Postnatal exposures only.** Pregnancy-window versions of the same families exist in
the release and are excluded deliberately: an exposure measured seven to ten years
before an assay is a different and much weaker question.

**Two things worth knowing before interpreting results.** The serum metabolome
features are indexed `metab_1` … `metab_177` with a chemical class and no compound
identifiers, because none exist upstream, so **loadings on that block are not
biologically interpretable**. And `bmi_cat` has only 9 children in its first level, so
it should not be used as a response.

**Cohort levels are the codes 1–6.** The release documents only "Cohort of inclusion
(1 to 6)" and nowhere maps those codes to the six countries, so the mapping has
deliberately not been guessed.

`Data/build_helix_data.R` produces the file from the prepared blocks and records every
step in `helix$provenance`.

---

## Licence, citation and acknowledgement

The HELIX release permits educational use: in the words of its own documentation,
*"These data may also be used for educational purposes."*

**A specific citation paragraph is required on any publication**, and it must be added
verbatim rather than paraphrased:

> This data were created as part of the ISGlobal Exposome data challenge 2021, presented in
> this publication (preprint: https://arxiv.org/abs/2202.01680 - under review in Env. Int.).
> The HELIX study [Vrijheid, Slama, et al. EHP 2014; Maitre et al. 2018 BMJ Open]
> represents a collaborative project across six established and ongoing longitudinal
> population-based birth cohort studies in six European countries (France, Greece,
> Lithuania, Norway, Spain, and the United Kingdom). The research leading to these results
> has received funding from the European Community's Seventh Framework Programme
> (FP7/2007-2013) under grant agreement no 308333 – the HELIX project and the
> H2020-EU.3.1.2. - Preventing Disease Programme under grant agreement no 874583 (ATHLETE
> project). The data used for the analyses described in this manuscript were obtained from:
> Figshare https://figshare.com/account/home#/projects/98813 (project number 98813 accesed
> on MM/DD/YYYY) and github
> https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021/.

### Acknowledgement

> **Draft — for review before publication.** This wording has not been agreed with
> ISGlobal and should be checked with Augusto Anguita before this repository is
> publicised.

The extract distributed in this repository is shared for workshop teaching with the
agreement of ISGlobal, with thanks to Augusto Anguita. We are grateful to the HELIX
consortium, and to the families who took part in the six cohorts, without whom none of
this material would exist.

---

## Building and rendering

The practical renders to HTML, which is what participants read, and the rendered file
is committed.

```sh
# pandoc is not always on PATH; a usable copy ships with Quarto
export RSTUDIO_PANDOC="/Applications/quarto/bin/tools"

Rscript -e 'rmarkdown::render("practical/exposome_analysis.Rmd", output_format = "html_document")'
```

### Showing or hiding the chunk output

The `01-options` chunk near the top of the `.Rmd` sets one flag, `show.results`,
which controls the output of every code chunk in the practical:

- `show.results <- FALSE` renders the code and nothing else. **This is the version
  participants read**: they run each line themselves and see the results in their own
  console. The committed HTML is built this way, and is about 1 MB.
- `show.results <- TRUE` also prints the results and draws the figures. Use it to
  check the material before a workshop. The HTML is then about 15 MB, because the
  twenty-five figures are embedded in it, so do not commit a render made this way.

Inline `` `r ...` `` results, such as the numbers quoted in the exercise answers, are
computed either way and always appear.

To rebuild the data set you need the prepared blocks from the scoping repository; see
the header of `Data/build_helix_data.R`.

```sh
Rscript Data/build_helix_data.R
```

### Provenance

The practical is adapted from the
[introduction to mixOmics workshop](https://github.com/mixOmics-org/mxoorg-workshop-intro-mixomics),
keeping its structure and much of its prose, with the methods changed from the
discriminant family to the regression family. The data preparation derives from the
Lê Cao Lab's scoping of the Exposome Data Challenge release, which records every
preparation decision and the reasoning behind it.

### Licence

Copyright © 2026 Kim-Anh Lê Cao. Licensed under AGPL-3.0-or-later.
See [LICENSE](LICENSE) for the full text.

The licence above covers the **practical material in this repository**. It does not
and cannot relicense the HELIX data, which remains subject to ISGlobal's terms as set
out above.
