# Vector heatmap panels for Figure 5D (plum hub genes) and Figure 7C (representative genes).
# Values are log2(TPM + 1) from data/Plum_hub_expression.csv and data/Representative_expression.csv;
# each row is z-standardised across the 36 libraries and truncated to [-2, 2] for display.
# Sampling labels follow the revised wording: 0 h (pre-injection) and 24 h post-challenge.
suppressPackageStartupMessages({library(ComplexHeatmap);library(circlize);library(grid)})

col_genotype_hm <- c(DD='#0072B2', GG='#D55E00', GD='#CC79A7')
col_age_hm <- c('1-year'='#FDB462', '2-year'='#B3B3B3')
time_levels <- c('0 h (pre-injection)', '24 h post-challenge')
col_time_hm <- setNames(c('#E8E8E8', '#333333'), time_levels)
col_z <- colorRamp2(c(-2, 0, 2), c('#2166AC', '#F7F7F7', '#B2182B'))

# Library names: lower case = 1-year, upper case = 2-year; dd = DD, ss = GG, sd = GD; _c_ = 0 h.
library_info <- function(ids) {
  data.frame(sample=ids,
             age=ifelse(grepl('^[a-z]', ids), '1-year', '2-year'),
             genotype=c(dd='DD', ss='GG', sd='GD')[tolower(substr(ids, 1, 2))],
             time=ifelse(grepl('_c_', ids, ignore.case=TRUE), time_levels[1], time_levels[2]),
             rep=as.integer(sub('.*_', '', ids)), stringsAsFactors=FALSE)
}

read_expr <- function(path) {
  x <- read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
  m <- as.matrix(x[, -(1:2)]); rownames(m) <- x$gene_id
  list(m=m, label=setNames(x$label, x$gene_id))
}

zscore <- function(m) { z <- t(scale(t(m))); z[z > 2] <- 2; z[z < -2] <- -2; z }

grid_mm <- unit(2.6, 'mm')
time_legend <- function(fs) list(title='Condition', title_gp=gpar(fontsize=fs, fontface='bold'),
                                 labels_gp=gpar(fontsize=fs - 0.5), legend_gp=gpar(col='#999999', lwd=.6),
                                 grid_height=grid_mm, grid_width=grid_mm)

top_anno <- function(info, order_names, fs, size) {
  HeatmapAnnotation(
    df=data.frame(Age=info$age, Genotype=factor(info$genotype, levels=c('DD','GG','GD')),
                  Condition=factor(info$time, levels=time_levels))[, order_names, drop=FALSE],
    col=list(Age=col_age_hm, Genotype=col_genotype_hm, Condition=col_time_hm),
    annotation_name_side='left', annotation_name_gp=gpar(fontsize=fs),
    simple_anno_size=size, border=TRUE, gp=gpar(col='#CCCCCC', lwd=.4),
    annotation_legend_param=list(
      Age=list(title_gp=gpar(fontsize=fs, fontface='bold'), labels_gp=gpar(fontsize=fs - 0.5), grid_height=grid_mm, grid_width=grid_mm),
      Genotype=list(title_gp=gpar(fontsize=fs, fontface='bold'), labels_gp=gpar(fontsize=fs - 0.5), grid_height=grid_mm, grid_width=grid_mm),
      Condition=time_legend(fs)))
}

z_legend <- function(fs, horizontal=FALSE) {
  l <- list(title='Z-score', title_gp=gpar(fontsize=fs, fontface='bold'), labels_gp=gpar(fontsize=fs - 0.5), at=c(-2, -1, 0, 1, 2))
  if (horizontal) c(l, list(direction='horizontal', legend_width=unit(18, 'mm'), grid_height=grid_mm))
  else c(l, list(grid_width=grid_mm))
}

# Figure 5D: ten highest-kME plum genes; kME and categories from Supplementary Table 3.
plum_hub_heatmap <- function(expr_csv, hub_csv, fs=7) {
  e <- read_expr(expr_csv)
  hub <- read.csv(hub_csv, skip=2, check.names=FALSE, stringsAsFactors=FALSE)
  hub <- hub[order(-hub$kME), ]
  stopifnot(nrow(hub) == 10, all(hub$Transcript_ID %in% rownames(e$m)))
  info <- library_info(colnames(e$m))
  ord <- order(info$age, match(info$genotype, c('DD','GG','GD')), match(info$time, time_levels), info$rep)
  info <- info[ord, ]; m <- e$m[hub$Transcript_ID, info$sample]
  z <- zscore(m)
  labels <- sub(' \\(kME=.*', '', hub$Hub_label)
  stopifnot(identical(unname(e$label[hub$Transcript_ID]), labels))
  right <- rowAnnotation(
    kME=anno_barplot(hub$kME, bar_width=.6, gp=gpar(fill='#888888', col=NA), ylim=c(0, 1),
                     width=unit(9, 'mm'), axis_param=list(side='top', at=c(0, .5, 1), labels=c('0', '0.5', '1'), labels_rot=0, gp=gpar(fontsize=fs - 1))),
    Category=anno_text(hub$Functional_category, gp=gpar(fontsize=fs - 0.5), location=0, just='left', width=unit(19, 'mm')),
    annotation_name_gp=gpar(fontsize=fs), annotation_name_rot=0, annotation_name_side='top', show_annotation_name=c(kME=TRUE, Category=FALSE))
  Heatmap(z, name='Z-score', col=col_z, cluster_rows=FALSE, cluster_columns=FALSE,
          row_labels=labels, row_names_side='left', row_names_gp=gpar(fontsize=fs, fontface='bold'),
          show_column_names=FALSE, column_split=factor(info$age, levels=c('1-year','2-year')), column_title=NULL,
          column_gap=unit(1, 'mm'), top_annotation=top_anno(info, c('Age','Genotype','Condition'), fs, unit(2.6, 'mm')),
          right_annotation=right, border_gp=gpar(col='#CCCCCC', lwd=.4),
          heatmap_legend_param=z_legend(fs, horizontal=TRUE))
}

# Figure 7C: four turquoise and four plum representative genes; kME from Supplementary Table 7.
representative_heatmap <- function(expr_csv, table7_csv, fs=7.5) {
  e <- read_expr(expr_csv)
  t7 <- read.csv(table7_csv, skip=2, check.names=FALSE, stringsAsFactors=FALSE)
  t7 <- t7[!is.na(t7$kME) & t7$gene_id != '', ]
  path_label <- c(CHST11='GAG biosynthesis', ODC1='Arg/Pro metabolism', GALNT='O-glycan biosynthesis',
                  'FUT1/2'='Glycosphingolipid', PRKCD='NOD-like/FcgammaR', NOTCH1='Notch signaling',
                  GSN='FcgammaR phagocytosis', CCNG1='p53 signaling')
  genes <- c('CHST11','ODC1','GALNT','FUT1/2','PRKCD','NOTCH1','GSN','CCNG1')
  t7 <- t7[match(genes, t7$gene_symbol), ]
  stopifnot(!anyNA(t7$gene_id), all(t7$gene_id %in% rownames(e$m)))
  info <- library_info(colnames(e$m))
  ord <- order(info$age, match(info$time, time_levels), match(info$genotype, c('DD','GG','GD')), info$rep)
  info <- info[ord, ]; m <- e$m[t7$gene_id, info$sample]
  z <- zscore(m)
  module <- factor(t7$module, levels=c('Turquoise','Plum'))
  row_lab <- sprintf('%s (%s)', t7$gene_symbol, path_label[t7$gene_symbol])
  left <- rowAnnotation(
    Module=module, col=list(Module=c(Turquoise='#40E0D0', Plum='#DDA0DD')), simple_anno_size=unit(2.6, 'mm'),
    show_annotation_name=FALSE, show_legend=TRUE,
    annotation_legend_param=list(Module=list(title_gp=gpar(fontsize=fs, fontface='bold'), labels_gp=gpar(fontsize=fs - 0.5), grid_height=grid_mm, grid_width=grid_mm)))
  right <- rowAnnotation(
    kME=anno_barplot(t7$kME, bar_width=.6, gp=gpar(fill='#888888', col=NA), ylim=c(0, 1), width=unit(12, 'mm'),
                     axis_param=list(side='top', at=c(0, .5, 1), labels=c('0', '0.5', '1'), labels_rot=0, gp=gpar(fontsize=fs - 1))),
    annotation_name_gp=gpar(fontsize=fs), annotation_name_rot=0, annotation_name_side='top')
  Heatmap(z, name='Z-score', col=col_z, cluster_rows=FALSE, cluster_columns=FALSE,
          row_labels=row_lab, row_names_side='right', row_names_gp=gpar(fontsize=fs),
          row_split=module, row_title_gp=gpar(fontsize=fs + 0.5, fontface='bold'), row_title_rot=90, row_gap=unit(1, 'mm'),
          show_column_names=FALSE, column_split=factor(info$age, levels=c('1-year','2-year')), column_title=NULL,
          column_gap=unit(1.2, 'mm'), top_annotation=top_anno(info, c('Age','Condition','Genotype'), fs, unit(3, 'mm')),
          left_annotation=left, right_annotation=right, border_gp=gpar(col='#CCCCCC', lwd=.4),
          heatmap_legend_param=z_legend(fs, horizontal=TRUE))
}

heatmap_grob <- function(ht, legend_side='right') {
  grid.grabExpr(draw(ht, heatmap_legend_side=legend_side, annotation_legend_side=legend_side, merge_legend=TRUE,
                     legend_gap=unit(2, 'mm'), padding=unit(c(2, 1, 7, 3), 'mm')))
}
