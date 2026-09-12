# Age-class-dependent survival heterosis in hybrid abalone

Analysis code and processed data accompanying the September 2026 revision of **Age-class-dependent survival heterosis in hybrid abalone challenged with Vibrio harveyi**.

The study compares *Haliotis gigantea* (GG), *H. discus hannai* (DD), and the female GG × male DD hybrid (GD) at one and two years of age. This release reproduces the reported analyses and figures from processed data. It is not a reconstruction of the raw sequencing or instrument-processing pipelines.

## Contents

- `data/`: current tables used by the analysis and figure scripts, exported without changing cell values. `manifest.json` identifies the source table, dimensions and SHA-256 digest of each export.
- `Analysis_Code/`: statistical recomputation and figure-generation scripts.
- `Analysis_Code/figure_inputs/`: the processed APA and histology inputs and only the image panels needed in the revised figures.

Earlier analysis versions, unused historical workbook sheets, reviewer correspondence and author annotations are excluded.

## Reproduce the results

Run the commands from the repository root. The Python scripts require NumPy, pandas, SciPy, openpyxl and Matplotlib. The R scripts require ggplot2, dplyr, tidyr, readr, readxl, patchwork, cowplot, magick and ragg. The R PDF device requires Cairo support.

```bash
python Analysis_Code/recompute_current.py
python Analysis_Code/recompute_plum_interaction.py
python Analysis_Code/redraw_apa_supplement.py

Rscript Analysis_Code/restore_current_phenotype_figures.R .
Rscript Analysis_Code/restore_molecular_figures.R .
Rscript Analysis_Code/restore_design_and_evidence.R Analysis_Code/figure_outputs
Rscript Analysis_Code/redraw_qpcr_profiles.R Source_Data_Revised.xlsx Analysis_Code/figure_outputs
```

The filename in the final command identifies the logical source workbook. When that workbook is absent, the reader uses the matching files in `data/`. Authors can also run the same scripts alongside the revision workbooks. Numeric recomputation results are written to `Analysis_Code/recomputed/`, and figures to `Analysis_Code/figure_outputs/`.

## Interpretation of the data

Survival contrasts use three independently supplied tanks per genotype–age group. Cardiac measurements are supplied as six group means; the cardiac resilience index is a descriptive principal-component summary, not a tank-level significance test. Low-dose cellular and histological time courses use different animals at successive times. Their trajectory integrals are descriptive group summaries.

RNA-seq expression is represented by 36 libraries, each prepared from a three-animal RNA pool. Module–trait comparisons use six genotype–age means, with unrestricted and age-stratified permutation references and correction across 246 comparisons. They do not pair transcriptomic and cardiac animals.

The RT-qPCR data contain 350 quality-filtered target records across 46 stored sample identifiers. Technical triplicates were averaged, and relative expression was calibrated to the gene-specific mean ΔCt of one-year-old GD. The identifier-to-time and RNA-pool assignments are unavailable. Accordingly, Figure 8 displays genotype–age RT-qPCR profiles beside the corresponding RNA-seq expression at 0 and 24 h, including both concordant and discordant patterns. This is a descriptive comparison, not time-matched validation. No individual RT-qPCR–cardiac association or inferred biological sample count is included. Additional primer matches are computational predictions, not empirical amplicon-specificity measurements.

The cardiac trace in Figure 2E is explicitly schematic and illustrates endpoint definitions; it is not an observed heart-rate recording. Histological images are preserved source panels. Lesion-positive proportions are reported descriptively because the underlying positive and total scoring counts are unavailable.

## Raw sequencing data

Raw sequencing data are deposited in the Genome Sequence Archive under **CRA042084** and are scheduled for public release after publication. They are not included in this repository. Reviewer access is handled separately by the corresponding author.

## Licence

The code retains the MIT licence used by the authors' earlier code repository. See `LICENSE`. This statement does not assign a separate open-data or image licence.
