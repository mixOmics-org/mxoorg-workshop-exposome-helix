## ---------------------------------------------------------------------------
## build_helix_data.R
##
## Builds Data/helix.RData, the single object the workshop practical loads.
##
## Input is the prepared HELIX blocks from the scoping repository,
## lecaolab-data-mxo, which records every preparation decision and why it was
## made. Nothing is re-derived here: this script selects, renames and reshapes
## what that repository already produced, so that the path from the ISGlobal
## release to the workshop extract is reproducible rather than manual.
##
## The scoping repository's prepared/ directory is gitignored, so it exists on
## disk but not in git. If the .RData files are missing, rebuild them:
##
##   cd ~/src/area/lecaolab/lecaolab-data-mxo/01-exposome-challenge
##   bash data/scripts/01_download_helix.sh
##   cd pre-processing
##   Rscript -e 'for (f in c("00-cohort","01-proteome","02-serum-metabolome",
##                           "03-urine-metabolome","04-transcriptome","05-exposome"))
##                 rmarkdown::render(paste0(f, ".Rmd"))'
##
## 00-cohort.Rmd must run first.
##
## Usage, from the repository root:
##   Rscript Data/build_helix_data.R
##   Rscript Data/build_helix_data.R /path/to/01-exposome-challenge/pre-processing/prepared
## ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

PREP <- if (length(args) >= 1) args[1] else path.expand(file.path(
  "~", "src", "area", "lecaolab", "lecaolab-data-mxo",
  "01-exposome-challenge", "pre-processing", "prepared"))

OUT <- file.path("Data", "helix.RData")

if (!dir.exists(PREP))
  stop("prepared/ not found at:\n  ", PREP,
       "\nPass its path as the first argument, or rebuild it (see the header).")

message("reading prepared blocks from: ", PREP)

src <- new.env()
for (f in c("helix_cohort.RData", "prot_prepared.RData", "serum_prepared.RData",
            "urine_prepared.RData", "gx_prepared.RData", "expo_prepared.RData")) {
  p <- file.path(PREP, f)
  if (!file.exists(p)) stop("missing input: ", p)
  load(p, envir = src)
}

cohort <- src$cohort
prot   <- src$prot
serum  <- src$serum
urine  <- src$urine
gx     <- src$gx
expo   <- src$expo

subjects <- cohort$subjects
meta     <- cohort$meta
codebook <- expo$codebook

stopifnot(length(subjects) == 896)

## ---------------------------------------------------------------------------
## 1. Omics blocks
##
## Four blocks, samples in rows, taken as prepared. Only the transcriptome is
## reduced further, and only for size -- see below.
## ---------------------------------------------------------------------------

## The scoping block is already 28,738 -> 5,000 (CallRate >= 80%, then
## variance). For the workshop we take the top 1,000 of those 5,000 by variance.
##
## The reason is size, not science. The full 5,000 makes this .RData about
## 34 MB; 1,000 brings it to roughly 7 MB, which matters when a room of people
## clones the repository at once. The 1,000 retain 55.6% of the 5,000-probe
## variance and every one of them has a gene symbol, so loadings stay
## interpretable. This is a teaching-material decision and is recorded as one.
v      <- apply(gx$data, 2, var)
keep   <- names(sort(v, decreasing = TRUE))[1:1000]
gx1000 <- gx$data[, keep, drop = FALSE]

var_retained <- sum(v[keep]) / sum(v)
message(sprintf("transcriptome: 5000 -> 1000 probes, %.1f%% of the variance retained",
                100 * var_retained))

omics <- list(
  proteome      = prot$data,
  serum         = serum$data,
  urine         = urine$data,
  transcriptome = gx1000)

## ---------------------------------------------------------------------------
## 2. Exposures: the nine POSTNATAL families
##
## Postnatal only. The pregnancy-window versions of the same families exist in
## the release and are excluded on purpose: an exposure measured seven to ten
## years before an assay is a different and much weaker question. Being
## measured at the same visit as the blood and urine is what gives these
## regressions a chance.
##
## Family membership comes from the codebook, never from parsing variable
## names. The friendly names on the left are the workshop's; the strings on the
## right are the literal codebook$family values.
## ---------------------------------------------------------------------------

FAMILIES <- c(
  metals            = "Metals",
  organochlorines   = "Organochlorines",
  phthalates        = "Phthalates",
  phenols           = "Phenols",
  pfas              = "Per- and polyfluoroalkyl substances (PFAS)",
  op_pesticides     = "Organophosphate pesticides",
  air_pollution     = "Air Pollution",
  built_environment = "Built environment",
  meteorological    = "Meteorological")

num <- expo$numeric
idx <- match(colnames(num), codebook$variable_name)
stopifnot(!anyNA(idx))

exposure <- lapply(FAMILIES, function(f) {
  sel <- codebook$family[idx] == f & codebook$period[idx] == "Postnatal"
  as.matrix(num[, sel, drop = FALSE])
})

## The expected widths, stated so a change upstream fails here rather than
## silently reshaping the practical.
expected_p <- c(metals = 9, organochlorines = 9, phthalates = 11, phenols = 7,
                pfas = 5, op_pesticides = 4, air_pollution = 12,
                built_environment = 15, meteorological = 9)
stopifnot(all(sapply(exposure, ncol)[names(expected_p)] == expected_p))

## The section 2 response must be present and in the family section 1 explores.
stopifnot("hs_sumPCBs5_cadj_Log2" %in% colnames(exposure$organochlorines))

## ---------------------------------------------------------------------------
## 3. Outcomes: categorical, all factors, all labelled
##
## Taken from cohort$meta rather than expo$phenotype. The two disagree:
## hs_asthma is numeric 0/1 in phenotype and already relabelled no/yes in meta,
## and meta is the only place the labels exist.
##
## Three of the five ship as bare numeric codes. Their labels are documented in
## codebook$description, so they are applied from the source rather than
## invented:
##   h_edumc_None  1: primary school, 2: secondary school, 3: university or higher
##   hs_bmi_c_cat  1: Thinness, 2: Normal, 3: Overweight, 4: Obese
##
## h_cohort is the exception. The codebook says only "Cohort of inclusion
## (1 to 6)" and nowhere in the release maps those codes to the six countries,
## so the levels are left as 1-6. Guessing the mapping would put invented
## provenance in a teaching dataset.
## ---------------------------------------------------------------------------

outcomes <- data.frame(
  sex     = meta$e3_sex_None,
  cohort  = meta$h_cohort,
  bmi_cat = factor(meta$hs_bmi_c_cat, levels = c(1, 2, 3, 4),
                   labels = c("thinness", "normal", "overweight", "obese")),
  asthma  = meta$hs_asthma,
  mat_edu = factor(meta$h_edumc_None, levels = c(1, 2, 3),
                   labels = c("primary", "secondary", "university")),
  row.names = subjects,
  stringsAsFactors = FALSE)

stopifnot(all(sapply(outcomes, is.factor)))

## ---------------------------------------------------------------------------
## 4. Covariates: numeric
## ---------------------------------------------------------------------------

covariates <- data.frame(
  child_age       = meta$hs_child_age_None,  # years, at examination
  zbmi_who        = meta$hs_zbmi_who,        # WHO BMI z-score
  birthweight     = meta$e3_bw,              # grams
  raven_score     = meta$hs_correct_raven,   # Raven's matrices, correct answers
  behaviour_total = meta$hs_Gen_Tot,         # CBCL total score
  row.names = subjects)

stopifnot(all(sapply(covariates, is.numeric)))

## ---------------------------------------------------------------------------
## 5. Annotation, each table keyed by its block's column names
##
## Keyed so that annotation$<block>[colnames(omics$<block>), ] resolves. The
## scoping repository had a bug where the urine table kept metab_N row names
## after the columns were renamed to compound names, which made every lookup by
## name return NA while positional access kept working. It is fixed upstream;
## the assertions below are here so it cannot come back through this script.
##
## The serum block carries a chemical class and no compound identifier, because
## none exists upstream. Its features are metab_1 ... metab_177 and its
## loadings are not biologically interpretable. That is a property of the
## release, recorded rather than worked around.
## ---------------------------------------------------------------------------

annotation <- list(
  proteome      = prot$fdata[, c("Prot_ID", "UniProtKB", "Gene_Symbol", "Gene_Name")],
  serum         = serum$fdata[, c("Class", "Class_2")],
  urine         = urine$fdata[, c("var", "CHEBI", "KEGG")],
  transcriptome = gx$fdata[keep, c("GeneSymbol_Affy", "EntrezeGeneID_Affy", "CallRate")],
  codebook      = codebook[codebook$variable_name %in% unlist(lapply(exposure, colnames)),
                           c("variable_name", "family", "subfamily", "description",
                             "transformation", "labelsshort")])

rownames(annotation$codebook) <- annotation$codebook$variable_name

## ---------------------------------------------------------------------------
## 6. Provenance
## ---------------------------------------------------------------------------

provenance <- list(
  source = paste("HELIX (Human Early-Life Exposome), released by ISGlobal for the",
                 "Exposome Data Challenge 2021. Six population-based birth cohorts",
                 "in France, Greece, Lithuania, Norway, Spain and the United Kingdom."),
  upstream = c(
    data  = "github.com/isglobal-exposomeHub/ExposomeDataChallenge2021",
    files = "github.com/isglobal-brge/brge_data_large (data/ExposomeDataChallenge2021, Git LFS)"),
  scoping = paste("Prepared blocks taken from lecaolab-data-mxo,",
                  "01-exposome-challenge/pre-processing/prepared/."),
  prepared = as.character(Sys.Date()),
  citation_required = paste(
    "Educational use is permitted. A specific citation paragraph is required",
    "verbatim on any publication; it is reproduced in README.md and in the",
    "practical, and must be added rather than paraphrased."),
  redistribution = paste(
    "Redistribution of this extract for the workshop was agreed with",
    "Augusto Anguita at ISGlobal."),
  steps = c(
    "n = 896: the intersection of all five blocks; the transcriptome binds",
    "all 896 samples in the same order in every element, asserted on build",
    "proteome 36, serum 177, urine 44: taken as prepared, no further filtering",
    "no transformation applied to any block: all four arrived transformed",
    "  (proteome and serum log, urine log2 with a 0.05 uM floor, transcriptome",
    "  log2 and per-gene centred), established from the data upstream",
    sprintf("transcriptome 5000 -> 1000 by variance, retaining %.1f%% of the 5000-probe variance",
            100 * var_retained),
    "  the transcriptome reduction is for repository size, not a scientific choice",
    "exposures: the 9 postnatal families only, 81 variables, selected by codebook",
    "  family and period; pregnancy-window versions excluded as a different question",
    "outcomes taken from the cohort metadata, where asthma is labelled no/yes;",
    "  maternal education and BMI category labelled from codebook descriptions",
    "cohort levels left as 1-6: the release does not map the codes to countries",
    "annotation keyed by each block's column names, asserted on build"))

## ---------------------------------------------------------------------------
## 7. Assemble, assert, save
## ---------------------------------------------------------------------------

helix <- list(
  omics      = omics,
  exposure   = exposure,
  outcomes   = outcomes,
  covariates = covariates,
  annotation = annotation,
  provenance = provenance)

## Every element on the same 896 children, in the same order. Asserted here
## rather than trusted from the preparation documents, because it is cheap and
## because everything downstream depends on it.
for (b in names(helix$omics))
  stopifnot(nrow(helix$omics[[b]]) == 896,
            identical(rownames(helix$omics[[b]]), subjects))
for (b in names(helix$exposure))
  stopifnot(nrow(helix$exposure[[b]]) == 896,
            identical(rownames(helix$exposure[[b]]), subjects))
stopifnot(identical(rownames(helix$outcomes), subjects),
          identical(rownames(helix$covariates), subjects))

## Every annotation table resolves against its block's column names.
for (b in c("proteome", "serum", "urine", "transcriptome")) {
  a <- helix$annotation[[b]]
  stopifnot(identical(rownames(a), colnames(helix$omics[[b]])),
            !any(is.na(rownames(a))))
}

## Every exposure variable has a codebook row with a description.
all_exp <- unlist(lapply(helix$exposure, colnames), use.names = FALSE)
stopifnot(all(all_exp %in% rownames(helix$annotation$codebook)),
          !any(is.na(helix$annotation$codebook[all_exp, "description"])))

## Numeric and complete.
stopifnot(all(sapply(helix$omics, function(x) is.numeric(x) && !anyNA(x))),
          all(sapply(helix$exposure, function(x) is.numeric(x) && !anyNA(x))))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
save(helix, file = OUT, compress = "xz")

message("\nwrote ", OUT, sprintf(" (%.1f MB)", file.size(OUT) / 1024^2))
message("\nomics:")
print(sapply(helix$omics, dim))
message("exposure:")
print(sapply(helix$exposure, ncol))
message("outcomes:")
print(summary(helix$outcomes))
