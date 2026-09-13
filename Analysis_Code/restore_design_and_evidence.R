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
         width=w, height=h, units='mm', device=grDevices::cairo_pdf, bg='white')
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised.svg')), p,
         width=w, height=h, units='mm', device=svglite::svglite, bg='white')
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised.png')), p,
         width=w, height=h, units='mm', dpi=600, device=ragg::agg_png, bg='white')
  ggsave(file.path(figdir, paste0('Figure_', n, '_revised.tiff')), p,
         width=w, height=h, units='mm', dpi=600, device=ragg::agg_tiff,
         compression='lzw', bg='white')
}

# Five animal allocations; cellular assays and histology share the same animals.
# RT-qPCR uses aliquots of the 24 h short-read RNA-seq pools.
if (length(args) < 2 || args[[2]] %in% c('1', 'all')) {
p <- base(108) + txt(50,104,'Experimental design and animal allocation',4.0,TRUE) +
  box(20,94,23,11,'#0072B2') + txt(20,96,'DD male',3.6,TRUE,'white') +
  txt(20,91.5,'H. discus hannai',3.0,FALSE,'white',TRUE) +
  txt(35,94,'×',4.6,TRUE) + box(50,94,23,11,'#D55E00') +
  txt(50,96,'GG female',3.6,TRUE,'white') + txt(50,91.5,'H. gigantea',3.0,FALSE,'white',TRUE) +
  link(62,94,68,94) + box(81,94,24,11,'#CC79A7') + txt(81,94,'GD hybrid (F₁)',3.6,TRUE) +
  txt(50,83,'Each genotype: one-year-old and two-year-old animals',3.3,TRUE) +
  link(20,80,20,76) + link(60,80,60,76) + link(90,80,90,76) +
  box(20,71,38,10,'#F4A259') + txt(20,71,'High-dose challenge\n10⁹ CFU/mL',3.3,TRUE) +
  box(60,71,38,10,'#61A6A6') + txt(60,71,'Low-dose challenge\n10⁶ CFU/mL',3.3,TRUE) +
  box(90,71,18,10,'#CEE5FB') + txt(90,71,'Pre-injection\nreference',3.15,TRUE)
centres <- c(10,30,50,70,90)
fills <- c('#FFF9E9','#FFF9E9','#F0FAF8','#F0FAF8','#F0F6FC')
heads <- c('Survival\nanalysis','Cardiac\nrecording','Flow cytometry\nand histology','Short-read\nRNA-seq','Iso-Seq\nreference')
animals <- c('30 challenged\n+ 10 controls','8 challenged\n+ 3 controls','30 incl. baseline\n+ 10 controls','18 animals\n9 at each time','6 animals')
units <- c(rep('per genotype–age\ngroup',4),'per genotype\n3 from each age')
assay <- c('Deaths recorded at\n0, 24, 48,\n72 and 96 h',
           'Heart rate recorded\ncontinuously\nafter injection',
           '0, 24, 48, 96\nand 168 h:\n5 animals/time\nhaemolymph,\nthen tissues',
           '0 and 24 h:\n3 RNA pools/time\n3 animals/pool\n36 libraries total',
           '0 h:\n1 library/genotype\n3 libraries total')
for (i in seq_along(centres)) {
  x <- centres[[i]]
  p <- p + link(x,66,x,63) + box(x,41.8,18,42.4,fills[[i]],'grey60',.4) +
    txt(x,59.3,heads[[i]],3.2,TRUE) +
    txt(x,52.3,'Animals',3.15,TRUE,'#365463') + txt(x,47.4,animals[[i]],3.2) +
    txt(x,41,units[[i]],3.1) +
    txt(x,34.3,'Sampling & assay',3.05,TRUE,'#365463') + txt(x,27.6,assay[[i]],3.1)
}
p <- p + link(70,20.6,70,17.5) + box(50,8.8,98,17.4,'#E8F2FB','grey60',.4) +
  txt(50,14.0,'RT-qPCR — 24 h RNA-pool aliquots',3.45,TRUE) +
  txt(50,9.9,'Source: 3 RNA pools per genotype–age group (18 pools total)',3.2) +
  txt(50,5.0,'Target genes: CHST11, ODC1, GALNT, FUT1/FUT2,\nPRKCD, NOTCH1, GSN and CCNG1',3.2)
savefig(p,1,190,185)
}

if (length(args) < 2 || args[[2]] %in% c('9', 'all')) {
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
}
cat('Selected design and evidence diagrams exported.\n')
