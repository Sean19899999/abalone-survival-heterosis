# Physiological and transcriptomic analysis of age-dependent survival heterosis

Analysis code and processed data accompanying the revisions (September and November 2026) of **Physiological and transcriptomic analysis of age-dependent survival heterosis in hybrid abalone**.

The study compares *Haliotis gigantea* (GG), *H. discus hannai* (DD), and the female GG × male DD hybrid (GD) at one and two years of age. This release reproduces the reported analyses and figures from processed data. It is not a reconstruction of the raw sequencing or instrument-processing pipelines.

## Contents

- `data/`: current tables used by the analysis and figure scripts, exported without changing cell values. `manifest.json` identifies the source table, dimensions and SHA-256 digest of each export.
- `Analysis_Code/`: statistical recomputation and figure-generation scripts.
- `Analysis_Code/figure_inputs/`: the processed APA and histology inputs and only the image panels needed in the revised figures.

Earlier analysis versions, unused historical workbook sheets, reviewer correspondence and author annotations are excluded.

## Reproduce the results

Run the commands from the repository root. The Python scripts require NumPy, pandas, SciPy, openpyxl, Matplotlib and Pillow. The R scripts require ggplot2, dplyr, tidyr, readr, readxl, patchwork, cowplot, magick, ragg and svglite, and the Bioconductor packages ComplexHeatmap and circlize. The R PDF device requires Cairo support; set `FIG_NO_PDF=1` to write only PNG and TIFF files.

```bash
python Analysis_Code/recompute_current.py
python Analysis_Code/recompute_plum_interaction.py
python Analysis_Code/redraw_apa_supplement.py

# Histology panels (Figure 3E-F, Supplementary Figure 6A-B): scale bars, then lesion arrows
python Analysis_Code/redraw_histology_scalebars.py Analysis_Code/figure_inputs/histology_scalebar_spec.json
python Analysis_Code/draw_histology_arrows.py Analysis_Code/figure_inputs/histology_scalebar_spec.json Analysis_Code/figure_inputs/histology_arrows.json Analysis_Code/figure_outputs/histology
cp Analysis_Code/figure_outputs/histology/histology_panels_EF_arrows.png Analysis_Code/figure_inputs/histology_panels_EF.png
python Analysis_Code/supplementary_fig7_lesion_trajectories.py . Analysis_Code/figure_outputs/Supplementary_Figure_7.png

Rscript Analysis_Code/restore_current_phenotype_figures.R .
Rscript Analysis_Code/restore_molecular_figures.R .
Rscript Analysis_Code/restore_design_and_evidence.R Analysis_Code/figure_outputs
Rscript Analysis_Code/redraw_qpcr_profiles.R Source_Data_Revised.xlsx Analysis_Code/figure_outputs
```

The filename in the final command identifies the logical source workbook. When that workbook is absent, the reader uses the matching files in `data/`. Authors can also run the same scripts alongside the revision workbooks. Numeric recomputation results are written to `Analysis_Code/recomputed/`, and figures to `Analysis_Code/figure_outputs/`.

## Interpretation of the data

Survival contrasts use three independently supplied tanks per genotype–age group. Cardiac measurements are supplied as six group means; Figure 2F directly displays their descriptive principal-component cardiac resilience scores, not a tank-level significance test. Both CBTB and CZTB are measured from injection: CZTB is the total time to cardiac arrest, not the interval after CBTB. Low-dose cellular and histological time courses use different animals at successive times. Their trajectory integrals are descriptive group summaries.

RNA-seq expression is represented by 36 libraries, each prepared from a three-animal RNA pool. Module–trait comparisons use six genotype–age means, with unrestricted and age-stratified permutation references and correction across 246 comparisons. They do not pair transcriptomic and cardiac animals.

RT-qPCR used aliquots of the 24 h digestive-gland RNA pools prepared for short-read RNA-seq. The source cohort supplied three three-animal pools per genotype–age group at 24 h. The table contains 350 quality-filtered target measurements under 46 stored aliquot identifiers; technical triplicates were averaged, and relative expression was calibrated to the gene-specific mean ΔCt of one-year-old GD. Figure 8 compares genotype–age RT-qPCR summaries with the 24 h RNA-seq profiles. Across eight genes, two age classes and three genotype pairs, 34 of 48 directions agree (15/16 for GD–GG, 14/16 for DD–GG and 5/16 for DD–GD). The script exports every comparison as `Figure_8_genotype_direction_comparison.csv`. Pool-level links for the individual aliquot identifiers are unavailable, so the RT-qPCR standard deviations describe variation among aliquot measurements within groups. Additional primer matches are computational predictions of sequence matches.

The cardiac trace in Figure 2E is explicitly schematic and illustrates endpoint definitions; it is not an observed heart-rate recording. Histological images are preserved source panels. In the November 2026 revision their microscope scale bars were redrawn for legibility: every bar represents 100 µm, and its length is converted from the calibration of the original microscope bar in the same image (`histology_scalebar_spec.json`); only the bottom-right label box changes. Lesion arrows (`histology_arrows.json`) mark examples of disrupted tissue architecture, vacuolation or haemocyte infiltration. The image scripts use the Noto Sans font path given in the spec. Lesion-positive proportions are the number of lesion-positive animals among the five sampled per genotype–age–time–tissue group; an animal was classified as lesion-positive when any of its four sections showed a lesion. They are reported descriptively.

## Raw sequencing data

Raw sequencing data are deposited in the Genome Sequence Archive under **CRA042084** and are scheduled for public release after publication. They are not included in this repository. Reviewer access is handled separately by the corresponding author.

## Licence

The code retains the MIT licence used by the authors' earlier code repository. See `LICENSE`. This statement does not assign a separate open-data or image licence.
