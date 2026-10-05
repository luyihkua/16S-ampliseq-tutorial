# 16S rRNA amplicon analysis in R (nf-core/ampliseq → phyloseq)

A tutorial for analysing nf-core/ampliseq v2.18.0 output in R with phyloseq,
microViz, microbiome, vegan, ANCOMBC and ggpicrust2. The example data are the
MiSeq SOP faecal samples from one mouse: Early (E01–E09, days 0–9) and Late
(L01–L10, days 141–150).

**[Open the tutorial](01_16S_tutorial.md)**

| Section | Content |
|---|---|
| 1 | nf-core/ampliseq processing, SILVA 138.2, DADA2 vs QIIME 2 taxonomy |
| 2–11 | DADA2 taxonomy: design, raw abundance, library size, rarefaction, alpha diversity, relative abundance, PCoA, PERMANOVA, ANCOM-BC2 with top-5 boxplots, summary |
| 12 | QIIME 2 taxonomy: relative abundance, PERMANOVA, ANCOM-BC2, comparison with DADA2 |
| 13 | PICRUSt2 / ggpicrust2 functional analysis (EC, KO, KEGG, MetaCyc) |

## Repository layout

```
01_16S_tutorial.Rmd    source (R Markdown)
01_16S_tutorial.md     rendered page shown on GitHub
data/                  sample design
results/               result tables (dada2/, qiime2/, picrust2/)
figures/               figures (dada2/, qiime2/, picrust2/)
```

## Re-render

```r
install.packages(c("rmarkdown", "knitr", "readr", "readxl", "dplyr", "tidyr"))
rmarkdown::render("01_16S_tutorial.Rmd", output_format = "github_document")
```
