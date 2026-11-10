# Second-revision readability fixes for panels that are preserved source rasters (Figures 4B, 5A, 5E).
# The plotted graphics (heatmap cells, dendrogram, network) are kept exactly as in the source panels;
# only their text is removed and redrawn as vector text at a legible size, with revised sampling labels.
# Pixel coordinates refer to the cached 600-dpi panels in Analysis_Code/figure_inputs.
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(magick);library(grid);library(ComplexHeatmap);library(circlize)})

fig_input <- function(root, name) image_read(file.path(root, 'Analysis_Code/figure_inputs', name))
white_out <- function(im, x, y, w, h) image_composite(im, image_blank(w, h, 'white'), offset=sprintf('+%d+%d', x, y))

# Figure 4B: k-means heatmap of the top 1,000 variable genes. Heatmap body and its annotation bars are kept;
# the legend column and the left annotation names are redrawn.
panel_4B <- function(root, fs=6.5) {
  src <- fig_input(root, 'panel_4_0_0.331_0.543_0.396.png')
  body <- image_crop(src, '1701x1589+0+0')
  body <- white_out(body, 0, 0, 110, 110)              # old panel tag
  body <- white_out(body, 300, 40, 300, 160)           # old annotation names (Age, Genotype, Condition)
  cbar <- image_crop(src, '243x30+1718+727')
  W <- 1701; H <- 1589; px_mm <- W / 78                # heatmap body displayed ~78 mm wide
  L <- 21 * px_mm; B <- 15 * px_mm                     # right legend column and bottom legend strip
  grid_mm <- unit(2.6, 'mm'); tg <- gpar(fontsize=fs, fontface='bold'); lg <- gpar(fontsize=fs - 0.5)
  mk <- function(labels, title, fill, ...) Legend(labels=labels, title=title, legend_gp=gpar(fill=fill), title_gp=tg, labels_gp=lg,
                                                  grid_height=grid_mm, grid_width=grid_mm, ...)
  clus <- mk(c('C1 (n=159)','C2 (n=84)','C3 (n=323)','C4 (n=395)','C5 (n=39)'), 'Gene cluster',
             c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854'))
  bottom <- packLegend(mk(c('1-year','2-year'), 'Age', c('#FDB462','#B3B3B3')),
                       mk(c('DD','GG','GD'), 'Genotype', c('#0072B2','#D55E00','#CC79A7')),
                       mk(c('0 h (pre-injection)','24 h post-challenge'), 'Condition', c('white','#333333'), border='grey40'),
                       direction='horizontal', column_gap=unit(5, 'mm'))
  clus_grob <- grid.grabExpr(draw(clus, x=unit(0, 'npc'), y=unit(1, 'npc'), just=c('left', 'top')))
  bottom_grob <- grid.grabExpr(draw(bottom, x=unit(0, 'npc'), y=unit(1, 'npc'), just=c('left', 'top')))
  x0 <- W + 1.5 * px_mm; cbw <- 17 * px_mm; cy1 <- H - 9 * px_mm; cy0 <- cy1 - 2.6 * px_mm
  tick_x <- x0 + (c(-3, -1, 0, 1, 3) + 3) / 6 * cbw
  ggplot() +
    annotation_raster(as.raster(body), xmin=0, xmax=W, ymin=0, ymax=H) +
    annotate('text', x=595, y=H - c(69, 116.5, 165), label=c('Age','Genotype','Condition'), hjust=1, size=(fs - 0.5) / .pt, fontface='bold') +
    annotate('text', x=x0, y=H - 4.5 * px_mm, label='Expression (z)', hjust=0, vjust=.5, size=fs / .pt, fontface='bold') +
    annotation_raster(as.raster(cbar), xmin=x0, xmax=x0 + cbw, ymin=cy0, ymax=cy1) +
    annotate('segment', x=tick_x, xend=tick_x, y=cy0, yend=cy0 - .8 * px_mm, linewidth=.3) +
    annotate('text', x=tick_x, y=cy0 - 1.2 * px_mm, label=c('-3','-1','0','1','3'), vjust=1, size=(fs - 0.5) / .pt) +
    annotation_custom(clus_grob, xmin=x0, xmax=W + L, ymin=cy0 - 32 * px_mm, ymax=cy0 - 6 * px_mm) +
    annotation_custom(bottom_grob, xmin=200, xmax=W + L, ymin=-B, ymax=-1.5 * px_mm) +
    coord_fixed(xlim=c(0, W + L), ylim=c(-B, H), expand=FALSE) + theme_void()
}

# Figure 5A: residual-network dendrogram and module colours. Axis ticks at 1.0 and 0.5 lie at rows 149.5 and 1128.5.
panel_5A <- function(root, fs=7) {
  src <- fig_input(root, 'panel_5_0_0_0.544_0.375.png')
  y_of <- function(px) 1 - (px - 149.5) / 1958
  den <- image_crop(src, '1820x1030+236+140'); bar <- image_crop(src, '1820x250+236+1170')
  pd <- ggplot(data.frame(x=c(0, 1), y=c(y_of(1170), y_of(140))), aes(x, y)) + geom_blank() +
    annotation_raster(as.raster(den), xmin=0, xmax=1, ymin=y_of(1170), ymax=y_of(140)) +
    scale_y_continuous(breaks=seq(.5, 1, .1), labels=sprintf('%.1f', seq(.5, 1, .1)), expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) + coord_cartesian(xlim=c(0, 1), ylim=c(y_of(1170), y_of(140))) +
    labs(x=NULL, y='Height') +
    theme_classic(base_size=fs) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank(),
                                        axis.text.y=element_text(size=fs, colour='black'), axis.title.y=element_text(size=fs + 0.5),
                                        axis.line.y=element_line(linewidth=.3), axis.ticks.y=element_line(linewidth=.3), plot.margin=margin(16, 2, 0, 2))
  pb <- ggplot(data.frame(x=c(0, 1), y=c(0, 1)), aes(x, y)) + geom_blank() +
    annotation_raster(as.raster(bar), xmin=0, xmax=1, ymin=0, ymax=1) +
    scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(breaks=.5, labels='Module', expand=c(0, 0)) +
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1)) + labs(x=NULL, y=NULL) +
    theme_void(base_size=fs) + theme(axis.text.y=element_text(size=fs + 0.5, hjust=1, margin=margin(r=2)), plot.margin=margin(0, 2, 2, 2))
  plot_grid(pd, pb, ncol=1, align='v', axis='lr', rel_heights=c(1030, 250 * 1.15))
}

# Figure 5E: plum hub network. Network graphics are kept and enlarged; both legends are redrawn below it.
panel_5E <- function(root, fs=7) {
  src <- fig_input(root, 'panel_5_0_0.699_0.598_0.301.png')
  net <- image_crop(src, '1672x1155+0+30')
  net <- white_out(net, 0, 0, 120, 140)                # old panel tag
  s <- 88 / (1672 / 600 * 25.4)                        # display scale relative to the 600-dpi source
  lw <- c(9, 16, 22) / 600 * 25.4 * s / 0.753          # measured key thickness (px) -> ggplot linewidth
  cats <- c('Immune signaling'='#E41A1C', 'Signaling'='#984EA3', 'Ubiquitin/ER'='#377EB8', 'Cell cycle'='#4DAF4A', 'Other'='#CCCCCC')
  lt <- theme_void() + theme(legend.position='bottom', legend.text=element_text(size=fs),
                             legend.margin=margin(0, 0, 0, 0), legend.key.size=unit(3, 'mm'),
                             legend.key.spacing.y=unit(.6, 'mm'), legend.title=element_text(size=fs, face='bold', margin=margin(b=1)))
  le <- get_legend(ggplot(data.frame(w=factor(c('0.1','0.2','0.3')), x=1:3, y=1), aes(x, y, linewidth=w, group=1)) +
                     geom_line(colour='grey60') + scale_linewidth_manual(values=lw, name='Edge weight') + lt +
                     theme(legend.key.width=unit(6, 'mm')) + guides(linewidth=guide_legend(nrow=1)))
  lc <- get_legend(ggplot(data.frame(cat=factor(names(cats), levels=names(cats)), x=1:5, y=1), aes(x, y, fill=cat)) +
                     geom_point(shape=21, size=2.6, stroke=.4, colour='grey30') + scale_fill_manual(values=cats, name='Category') + lt +
                     guides(fill=guide_legend(nrow=2, byrow=TRUE, title.position='top')))
  net_h <- 1155 / 600 * 25.4 * s
  plot_grid(ggdraw() + draw_image(net), ggdraw(le), ggdraw(lc), ncol=1, rel_heights=c(net_h, 5.5, 12))
}
