args <- commandArgs(trailingOnly = TRUE)
figdir <- args[[1]]
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({library(ggplot2); library(grid)})
base <- function(ymax = 100) ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, ymax), expand = FALSE, clip = 'off') +
  theme_void() + theme(text = element_text(family = 'sans'), plot.margin = margin(4, 4, 4, 4))
box <- function(x, y, w, h, fill, border = 'grey35', lw = .55)
  annotate('rect', xmin = x-w/2, xmax = x+w/2, ymin = y-h/2, ymax = y+h/2,
           fill = fill, colour = border, linewidth = lw)
txt <- function(x, y, label, size = 3, bold = FALSE, colour = 'black', italic = FALSE)
  annotate('text', x = x, y = y, label = label, size = size,
           fontface = if (italic) 'italic' else if (bold) 'bold' else 'plain',
           colour = colour, lineheight = 1.08)
link <- function(x, y, xend, yend, direction = TRUE, linetype = 'solid')
  annotate('segment', x=x, y=y, xend=xend, yend=yend, colour='grey40',
           linewidth=.45, linetype=linetype,
           arrow=if(direction) arrow(length=unit(1.6,'mm')) else NULL)
savefig <- function(p, n, w, h) {
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised_preview.png')), p,
         width=w, height=h, units='mm', dpi=180, device=ragg::agg_png, bg='white')
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised.pdf')), p,
         width=w, height=h, units='mm', device=cairo_pdf, bg='white')
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised.png')), p,
         width=w, height=h, units='mm', dpi=600, device=ragg::agg_png, bg='white')
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised.tiff')), p,
         width=w, height=h, units='mm', dpi=600, device=ragg::agg_tiff,
         compression='lzw', bg='white')
}

# Preserve the established DD/GG/GD and challenge colours, sans font and box style.
p <- base(130) + txt(50,127,'Experimental design and animal allocation',4.1,TRUE) +
  box(20,116,24,14,'#0072B2') + txt(20,119,'DD male',4.1,TRUE,'white') +
  txt(20,113,'H. discus hannai',2.8,FALSE,'white',TRUE) +
  txt(35,116,'×',5.3,TRUE) + box(50,116,24,14,'#D55E00') +
  txt(50,119,'GG female',4.1,TRUE,'white') + txt(50,113,'H. gigantea',2.8,FALSE,'white',TRUE) +
  link(62,116,68,116) + box(82,116,28,14,'#CC79A7') +
  txt(82,118,'GD hybrid',4.1,TRUE) + txt(82,112.5,'GG female × DD male',2.7) +
  box(50,101,90,10,'#F2F2F2') +
  txt(50,101,'DD, GG and GD × one- and two-year-old animals\nSix 100-L PVC tanks; approximately 100 animals per tank\nApproximately 35 days of acclimation; 25 °C after initial warming',3.05,TRUE) +
  link(35,96,25,92) + link(65,96,75,92) +
  box(25,85,44,14,'#F4A259') + txt(25,85,'High-dose challenge\n10⁹ CFU/mL; 50 μL',3.6,TRUE) +
  box(75,85,44,14,'#61A6A6') + txt(75,85,'Low-dose challenge\n10⁶ CFU/mL; 50 μL',3.6,TRUE) +
  link(19,78,13.5,74) + link(31,78,36.5,74) +
  link(67,78,62,74) + link(83,78,88,74) +
  box(13.5,63,21,22,'#FFFFCC') +
  txt(13.5,65,'Survival',3.3,TRUE) + txt(13.5,59,'30 animals/group\n3 tanks × 10\n0–96 h',2.85) +
  box(36.5,63,21,22,'#FFFFCC') +
  txt(36.5,68,'Cardiac profiles',3.2,TRUE) +
  txt(36.5,60,'8 separate\nanimals/group\nContinuous recording\nCBTB / CZTB',2.7) +
  box(62,63,26,22,'#FFFFCC') +
  txt(62,69,'Flow cytometry\nand histology',3.15,TRUE) +
  txt(62,60,'30 animals/group\nincluding 0 h\n5 animals per time\nBlood then tissues',2.85) +
  box(88,63,20,22,'#CEE5FB') +
  txt(88,69,'RNA-seq',3.25,TRUE) +
  txt(88,60,'18 separate\nanimals/group\n0 / 24 h: 9 each\n3 pools × 3 per time',2.65) +
  box(25,43,44,12,'#F4F4F4') +
  txt(25,43,'Survival controls\n10 sterile-seawater-injected animals/group',2.8) +
  box(62,43,26,12,'#F4F4F4') +
  txt(62,43,'Histology references\n10 saline-injected animals/group\nMatched sampling times',2.5) +
  link(88,52,88,48) + box(88,40,20,16,'#CEE5FB') +
  txt(88,43,'RT-qPCR',3.1,TRUE) + txt(88,37,'Same RNA\n0 h and 24 h',2.9) +
  txt(38,30.5,'Cellular and histological sampling: 0, 24, 48, 96 and 168 h\nFour sections per animal; representative H&E images at 168 h',2.7) +
  box(50,18,90,15,'#CEE5FB') +
  txt(50,22,'Separate baseline cohort for Iso-Seq reference libraries',3.4,TRUE) +
  txt(50,16,'3 genotype-specific libraries; 6 separate pre-injection animals per library\n3 one-year + 3 two-year animals of the same genotype; independent RNA extraction',2.8) +
  txt(50,5,'Animal counts are per genotype–age group unless indicated.\nRT-qPCR reuses RNA-seq material and requires no additional animals.',2.65,FALSE,'grey30')
savefig(p,1,180,220)

turq <- '#1FB5AF'; plum <- '#BB50D2'
p <- base() + box(50,89,43,16,'#879FAA','#536F7B',1) +
  annotate('text',x=50,y=92,label='bolditalic(V.~harveyi)~bold(challenge)',parse=TRUE,size=5,colour='white') +
  txt(50,85,'Low-dose transcriptomic profiles\n0 h versus 24 h',3.5,FALSE,'white') +
  box(25,65,47,20,'#58DCD1',turq,1) + txt(25,70,'Turquoise module',4.8,TRUE) +
  txt(25,61,'1,163 genes; 259 OD-high genes at two years\nOD-high enrichment: 4.60-fold; FDR = 3.7 × 10⁻¹¹³',3.1) +
  box(75,65,47,20,'#E0AEE3',plum,1) + txt(75,70,'Plum module',4.8,TRUE) +
  txt(75,61,'42 genes; immune- and stress-related annotations\nNominal cardiac associations across six groups',3.1) +
  box(25,42,47,20,'#EEFFFF',turq,.65) +
  txt(25,47,'Candidate tissue-maintenance functions',3.4,TRUE,'#007477') +
  txt(25,39,'Glycan biosynthesis / extracellular matrix metabolism\nCHST11, FUT1/FUT2, ODC1, GALNT\nCoordinated biosynthetic expression',3.1,FALSE,'#007477') +
  box(75,42,47,20,'#FFF1FF',plum,.65) +
  txt(75,47,'Candidate cellular stress regulation',3.4,TRUE,'#890092') +
  txt(75,39,'Signalling / autophagy / cell-state functions\nPRKCD, NOTCH1, GSN, CCNG1\nCZTB: rho = −0.886; FDR = 1.00',3.1,FALSE,'#890092') +
  box(50,20,64,12,'#276FC0','#124EAD',1) +
  txt(50,22,'Observed high-dose survival heterosis',4.2,TRUE,'white') +
  txt(50,17,'Two-year GD: +40.0 percentage points over the mid-parent value',3.2,FALSE,'white') +
  txt(50,7,'Separate high-dose phenotype and low-dose transcriptome cohorts.\nDTU: 62 response-switching genes at two years. APA: gene-specific hybrid deviations.',2.9,FALSE,'grey30') +
  annotate('segment',x=c(40,60),xend=c(25,75),y=c(81,81),yend=c(76,76),
           colour=c(turq,plum),linewidth=.8,linetype='dashed')
savefig(p,9,180,115)
cat('Updated design and evidence diagrams exported. No RT-qPCR validation data were inferred.\n')
