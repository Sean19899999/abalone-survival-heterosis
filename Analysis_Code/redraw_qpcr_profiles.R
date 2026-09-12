args <- commandArgs(trailingOnly=TRUE)
wb <- args[[1]]
out <- args[[2]]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readxl);library(cowplot)})
source(file.path(dirname(wb),'Analysis_Code/sheet_io.R'))
pal <- c(DD='#0072B2', GG='#D55E00', GD='#CC79A7')
gorder <- c('PRKCD','NOTCH1','GSN','CCNG1','CHST11','ODC1','GALNT','FUT1/2')
groups <- c('DD-1-year','DD-2-year','GD-1-year','GD-2-year','GG-1-year','GG-2-year')
glab <- c('DD\n1y','DD\n2y','GD\n1y','GD\n2y','GG\n1y','GG\n2y')
theme_set(theme_bw(8)+theme(text=element_text(family='sans'), panel.grid=element_blank(),
 axis.text=element_text(size=6), axis.title=element_text(size=7),
 plot.title=element_text(size=8,hjust=.5), plot.margin=margin(5,3,3,3),
 legend.title=element_blank(), legend.text=element_text(size=8)))
qp <- read_sheet(wb,sheet='Current_RTqPCR_Observations')
stopifnot(nrow(qp)==350, !any(grepl('CZTB',names(qp))),length(unique(qp$stored_sample_id))==46)
qp <- qp %>% mutate(group=factor(group,levels=groups))
qplots <- lapply(gorder,function(g){
  z <- filter(qp,gene==g)
  sm <- z %>% group_by(group,genotype) %>% summarise(m=mean(relative_expression),s=sd(relative_expression),.groups='drop')
  ggplot(z,aes(group,relative_expression,fill=genotype))+
    geom_col(data=sm,aes(y=m),width=.64,colour='grey35',linewidth=.25,alpha=.55)+
    geom_errorbar(data=sm,aes(y=m,ymin=m-s,ymax=m+s),width=.18,linewidth=.35)+
    geom_point(position=position_jitter(width=.10,height=0,seed=123),size=1.05,shape=21,colour='grey30',stroke=.2)+
    scale_fill_manual(values=pal,guide='none')+scale_x_discrete(labels=glab)+
    scale_y_continuous(expand=expansion(mult=c(0,.10)))+
    labs(title=paste0(g,' (',nrow(z),' records)'),x=NULL,y='Relative expression')
})
meta <- read_sheet(wb,sheet='WGCNA_library_eigengenes') %>% select(sample_id,age,genotype,condition) %>%
 mutate(matrix_sample_id=sub('^GD_','SD_',sub('^gd_','sd_',sub('^GG_','SS_',sub('^gg_','ss_',sample_id)))),
        time=factor(ifelse(condition=='Control','0 h','24 h'),levels=c('0 h','24 h')),
        group=factor(paste(genotype,age,sep='-'),levels=groups))
rep <- read_sheet(wb,sheet='Representative_expression')
stopifnot(nrow(rep)==8,nrow(meta)==36,setequal(rep$label,gorder),setequal(setdiff(names(rep),c('gene_id','label')),meta$matrix_sample_id))
expr <- rep %>% pivot_longer(-c(gene_id,label),names_to='matrix_sample_id',values_to='log2_tpm_plus_1') %>%
  left_join(meta,by='matrix_sample_id') %>% mutate(xpos=as.numeric(group)+ifelse(time=='0 h',-.12,.12))
stopifnot(nrow(expr)==288,!anyNA(expr$group))
rplots <- lapply(gorder,function(g){
  z <- filter(expr,label==g)
  sm <- z %>% group_by(group,genotype,time,xpos) %>% summarise(m=mean(log2_tpm_plus_1),s=sd(log2_tpm_plus_1),n=n(),.groups='drop')
  stopifnot(all(sm$n==3))
  ggplot(z,aes(xpos,log2_tpm_plus_1,colour=genotype,shape=time))+
    geom_point(position=position_jitter(width=.04,height=0,seed=123),size=1.25,alpha=.65,stroke=.45)+
    geom_errorbar(data=sm,aes(y=m,ymin=m-s,ymax=m+s),width=.13,linewidth=.35,alpha=.9)+
    geom_point(data=sm,aes(y=m),size=1.8,stroke=.6)+
    scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c('0 h'=1,'24 h'=16))+
    scale_x_continuous(breaks=1:6,labels=glab,limits=c(.5,6.5))+
    scale_y_continuous(expand=expansion(mult=c(.02,.10)))+
    labs(title=g,x=NULL,y=expression(log[2](TPM+1)))+theme(legend.position='none')
})
legendplot <- rplots[[1]]+theme(legend.position='bottom')+guides(shape=guide_legend(override.aes=list(colour='black',alpha=1,size=2.4)))
leg <- get_legend(legendplot)
heading <- function(tag,label) ggdraw()+draw_label(tag,x=0,y=.5,hjust=0,size=14)+draw_label(label,x=.065,y=.5,hjust=0,size=10)
fig <- plot_grid(heading('A','RT-qPCR expression across genotype-age groups'),
                plot_grid(plotlist=qplots,ncol=4),
                heading('B','RNA-seq expression at 0 h and 24 h'),
                plot_grid(plotlist=rplots,ncol=4),leg,
                ncol=1,rel_heights=c(.16,2,.16,2,.16))
for(ext in c('png','pdf','tiff')){
  dev <- switch(ext,png=ragg::agg_png,pdf=cairo_pdf,tiff=ragg::agg_tiff)
  if(ext=='tiff')ggsave(file.path(out,paste0('Figure_8_revised.',ext)),fig,width=180,height=240,units='mm',dpi=600,device=dev,compression='lzw',bg='white')
  else ggsave(file.path(out,paste0('Figure_8_revised.',ext)),fig,width=180,height=240,units='mm',dpi=600,device=dev,bg='white')
}
ggsave(file.path(out,'Figure_8_revised_preview.png'),fig,width=180,height=240,units='mm',dpi=180,device=ragg::agg_png,bg='white')
write.csv(expr %>% select(gene_id,label,sample_id,age,genotype,time,group,log2_tpm_plus_1),file.path(out,'Figure_8_RNAseq_source.csv'),row.names=FALSE)
cat('Figure 8 exported: 350 unchanged qPCR records and 288 RNA-seq feature-library values. No qPCR time assignments or correlation tests.\n')
