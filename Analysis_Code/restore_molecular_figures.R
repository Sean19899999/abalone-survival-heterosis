args<-commandArgs(trailingOnly=TRUE);root<-args[[1]];dest<-root;fd<-file.path(dest,'Analysis_Code/figure_outputs');wb<-file.path(dest,'Source_Data_Revised.xlsx')
dir.create(file.path(root,'Analysis_Code/figure_outputs'),showWarnings=FALSE,recursive=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr);library(readxl);library(cowplot);library(patchwork);library(magick)})
pal<-c(DD='#0072B2',GG='#D55E00',GD='#CC79A7');ages<-c('1-year','2-year')
theme_set(theme_bw(8)+theme(text=element_text(family='sans'),axis.text=element_text(size=7),axis.title=element_text(size=8),strip.background=element_blank(),strip.text=element_text(face='bold',size=9),panel.grid.minor=element_blank(),panel.grid.major=element_line(colour='#EBEBEB',linewidth=.3),plot.title=element_text(size=9,face='bold'),legend.title=element_text(size=7,face='bold'),legend.text=element_text(size=7),plot.margin=margin(8,4,3,6)))
source(file.path(root,'Analysis_Code/sheet_io.R'))
tab<-function(n) read_sheet(wb,sheet=n)
savefig<-function(p,n,w=180,h=200){
 for(ext in c('png','pdf','tiff')){
  dev<-switch(ext,png=ragg::agg_png,pdf=cairo_pdf,tiff=ragg::agg_tiff)
  if(ext=='tiff')ggsave(file.path(fd,paste0('Figure_',n,'_revised.',ext)),p,width=w,height=h,units='mm',dpi=600,device=dev,compression='lzw',bg='white')
  else ggsave(file.path(fd,paste0('Figure_',n,'_revised.',ext)),p,width=w,height=h,units='mm',dpi=600,device=dev,bg='white')
 }
 ggsave(file.path(fd,paste0('Figure_',n,'_revised_preview.png')),p,width=w,height=h,units='mm',dpi=180,device=ragg::agg_png,bg='white')
}
crop_panel<-function(n,x,y,w,h){
 cached<-file.path(root,'Analysis_Code/figure_inputs',paste0('panel_',paste(c(n,x,y,w,h),collapse='_'),'.png'))
 if(file.exists(cached)) return(ggdraw()+draw_image(as.raster(image_read(cached)),x=0,y=0,width=1,height=1))
 im<-image_read(file.path(root,'Analysis_Code/figure_inputs',paste0('image',n,'.tif')));i<-image_info(im)
 im<-image_crop(im,geometry_area(round(i$width*w),round(i$height*h),round(i$width*x),round(i$height*y)))
 ggdraw()+draw_image(as.raster(im),x=0,y=0,width=1,height=1)
}
tag<-function(p,label) ggdraw(p)+draw_label(label,x=0,y=1,hjust=0,vjust=1,size=15,fontface='plain')
pca<-tab('PCA_recomputed')%>%mutate(genotype=factor(genotype,levels=c('DD','GG','GD')),age=factor(age,levels=1:2,labels=ages),condition=factor(condition,levels=c('0 h','24 h')))
pA<-ggplot(pca,aes(PC1,PC2,colour=genotype,size=age,shape=condition))+geom_point(stroke=.6)+scale_colour_manual(values=pal,name='Genotype')+scale_size_manual(values=c(1.7,2.7),name='Age')+scale_shape_manual(values=c(1,16),name='Condition')+labs(x='PC1 (34.4% variance)',y='PC2 (9.9% variance)')+theme(legend.position='right')
modes<-tab('OD_mode_summary')%>%mutate(mode=trimws(mode),mode=recode(mode,'Additive'='Intermediate/residual','True Over-dominance (High)'='OD-high','True Over-dominance (Low)'='OD-low'),mode=factor(mode,levels=c('OD-low','Low-Parent Dominance','Intermediate/residual','High-Parent Dominance','OD-high')))
modepal<-c('Intermediate/residual'='#377EB8','High-Parent Dominance'='#FF7F00','Low-Parent Dominance'='#4DAF4A','OD-high'='#E41A1C','OD-low'='#984EA3')
pC<-ggplot(modes,aes(age,percent,fill=mode))+geom_col(width=.62,colour='grey40',linewidth=.2)+geom_text(aes(label=ifelse(percent>4,sprintf('%.1f%%',percent),'')),position=position_stack(vjust=.5),size=2.4,colour='white',fontface='bold')+scale_fill_manual(values=modepal,name='Response category')+labs(x='Age class',y='Proportion of genes (%)')+theme(legend.position='bottom',legend.key.height=unit(2.5,'mm'),legend.key.width=unit(3,'mm'))+guides(fill=guide_legend(ncol=1))
od<-tab('OD_classification')%>%mutate(mode=trimws(mode),display=case_when(mode=='True Over-dominance (High)'~'OD-high',mode=='True Over-dominance (Low)'~'OD-low',TRUE~'Other'))
scatter<-function(ageval)ggplot(filter(od,age==ageval),aes(MP_response,H_response,colour=display))+geom_point(size=.35,alpha=.3)+geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey40')+scale_colour_manual(values=c('Other'='grey75','OD-high'='#E41A1C','OD-low'='#984EA3'),guide='none')+coord_cartesian(xlim=c(-15,15),ylim=c(-15,15))+labs(title=ageval,x='Mid-parent response log₂(24 h / 0 h)',y='Hybrid response log₂(24 h / 0 h)')+theme(plot.title=element_text(hjust=.5))
# 原稿 B 热图和 D 重叠图只整体裁取，数值、配色及注释不变。
pB<-crop_panel(4,0,.331,.543,.396)
pD<-crop_panel(4,0,.738,.334,.254)
go_panel<-function(sheet,ageval){
 g<-read_sheet(wb,sheet=sheet,skip=2)%>%arrange(padj)%>%slice_head(n=6)
 g$short<-vapply(g$Description,function(s){parts<-strwrap(s,24);if(length(parts)>2)parts<-c(parts[1],paste0(substr(parts[2],1,21),'...'));parts<-ifelse(nchar(parts)>30,paste0(substr(parts,1,27),'...'),parts);paste(parts,collapse='\n')},character(1))
 g$short<-factor(g$short,levels=rev(unique(g$short)))
 ggplot(g,aes(-log10(padj),short,size=target_count,fill=fold_enrichment))+geom_point(shape=21,colour='white',stroke=.2)+scale_fill_gradient(low='#FCA082',high='#67001F',name='Fold enrichment')+scale_size_area(max_size=4,name='Gene count')+labs(title=ageval,x='−log₁₀(FDR)',y=NULL)+theme(plot.title=element_text(size=8,hjust=.5),axis.text.y=element_text(size=6),axis.text.x=element_text(size=6),legend.position='right',legend.key.height=unit(2,'mm'),legend.key.width=unit(2,'mm'),legend.text=element_text(size=5),legend.title=element_text(size=5))
}
pF<-plot_grid(go_panel('Supp_Fig_1B_GO_1yr','1-year-specific OD-high'),go_panel('Supp_Fig_1B_GO_2yr','2-year-specific OD-high'),ncol=1)
pE<-plot_grid(scatter('1-year'),scatter('2-year'),ncol=2)
fig4<-plot_grid(plot_grid(tag(pA,'A'),tag(pC,'C'),ncol=2),plot_grid(pB,tag(pF,'F'),ncol=2,rel_widths=c(1.1,.9)),plot_grid(pD,tag(pE,'E'),ncol=2,rel_widths=c(.62,1.38)),ncol=1,rel_heights=c(1.05,1.45,.8))
savefig(fig4,4,180,224)

# 表达值已是 log2(TPM+1)。只进行每行 z 标准化；样本顺序沿用原年龄/处理/基因型分组。
eig<-tab('WGCNA_library_eigengenes')%>%mutate(age=factor(age,levels=ages),genotype=factor(genotype,levels=c('DD','GG','GD')),condition=factor(condition,levels=c('Control','Infected')))%>%arrange(age,condition,genotype,sample_id)
sample_order<-eig$sample_id
# 已有分析命名 GG/GD 对应表达矩阵旧名 SS/SD；不改变 0/24 h 条件。
sample_order<-sub('^gg_','ss_',sample_order);sample_order<-sub('^gd_','sd_',sample_order)
sample_order<-sub('^GG_','SS_',sample_order);sample_order<-sub('^GD_','SD_',sample_order)
stopifnot(length(sample_order)==36,!anyDuplicated(sample_order))
heat<-function(df,label_order=NULL){
 if(is.null(label_order))label_order<-df$label
 z<-as.matrix(df[,sample_order]);z<-t(scale(t(z)));z[z>2]<-2;z[z< -2]<- -2;colnames(z)<-sample_order
 z<-as.data.frame(z);z$label<-df$label;z<-pivot_longer(z,-label,names_to='sample',values_to='z')%>%mutate(sample=factor(sample,levels=sample_order),label=factor(label,levels=rev(label_order)))
 p<-ggplot(z,aes(sample,label,fill=z))+geom_tile()+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-2,2),name='Z-score')+labs(x=NULL,y=NULL)+scale_x_discrete(labels=function(s)ifelse(seq_along(s)%%3==2,s,''))+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=90,hjust=1,size=5.5),axis.text.y=element_text(face='bold',size=7),legend.key.height=unit(4,'mm'),legend.key.width=unit(2.5,'mm'))
 p
}
tests<-tab('WGCNA_1e9_all_tests');sel<-union(unique(tests$module[tests$p_exact_720<.1]),'turquoise');assoc<-filter(tests,module%in%sel)%>%mutate(module=factor(module,levels=rev(sort(sel))),trait=factor(trait,levels=c('Survival_AUC','CBTB','CZTB')))
p5B<-ggplot(assoc,aes(trait,module,fill=rho))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=sprintf('%.2f',rho)),size=3)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',limits=c(-1,1),name='Spearman\nrho')+scale_x_discrete(labels=c('Survival AUC','CBTB','CZTB'))+labs(x=NULL,y=NULL,subtitle='Six group means; all FDR ≥ 0.05')+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=25,hjust=1),plot.subtitle=element_text(size=7))
p5C<-ggplot(eig,aes(genotype,plum,fill=condition))+geom_boxplot(position=position_dodge(width=.7),width=.55,outlier.shape=NA,linewidth=.35)+geom_point(position=position_jitterdodge(jitter.width=.08,dodge.width=.7,seed=10),colour='grey40',size=1.2,alpha=.7)+facet_wrap(~age)+scale_fill_manual(values=c(Control='white',Infected='grey35'),labels=c('0 h','24 h'),name='Condition')+labs(x='Genotype',y='Plum eigengene')+theme(legend.position='bottom')
hub<-tab('Plum_hub_expression');p5D<-crop_panel(5,.502,.399,.498,.286)
enrich1<-tab('OD_module_enrichment_1yr')%>%arrange(desc(fold_enrichment))%>%slice_head(n=10)%>%mutate(module=factor(module,levels=rev(module)))
p5F<-ggplot(enrich1,aes(fold_enrichment,module,fill=as.character(module)))+geom_col(width=.7,colour='grey60',linewidth=.2)+scale_fill_identity()+labs(x='OD-high fold enrichment (1-year)',y=NULL)
p5A<-crop_panel(5,0,0,.544,.375);p5E<-crop_panel(5,0,.699,.598,.301)
fig5<-plot_grid(p5A,tag(p5B,'B'),tag(p5C,'C'),p5D,p5E,tag(p5F,'F'),ncol=2,rel_heights=c(1,1,1))
savefig(fig5,5,180,204)

rs<-data.frame(pattern=c('Hybrid-specific shift','Direction reversal','Over-response','Under-response'),count=c(35,23,3,1),colour=c('#1B9E77','#377EB8','#E97848','#999999'))
rs$pattern<-factor(rs$pattern,levels=rs$pattern)
q6A<-ggplot(rs,aes(pattern,count,fill=colour))+geom_col(width=.65,colour='grey40',linewidth=.3)+geom_text(aes(label=sprintf('%d\n(%.1f%%)',count,100*count/62)),vjust=-.2,size=2.8)+scale_fill_identity()+scale_y_continuous(limits=c(0,42))+scale_x_discrete(labels=c('Hybrid-specific\nshift','Direction\nreversal','Over-\nresponse','Under-\nresponse'))+labs(x='Response pattern',y='Number of RS genes')
ix<-data.frame(module=c('turquoise','plum','darkviolet','lightpink4'),count=c(0,0,2,2));ix$module<-factor(ix$module,levels=ix$module)
q6B<-ggplot(ix,aes(module,count,fill=as.character(module)))+geom_col(width=.65,colour='grey40',linewidth=.3)+geom_text(aes(label=count),vjust=-.5,size=3)+scale_fill_identity()+scale_y_continuous(limits=c(0,3),breaks=0:3)+labs(x='WGCNA module',y='Mapped RS features')+theme(axis.text.x=element_text(angle=25,hjust=1))
iso<-tab('DTU_UBOX5_libraries')%>%mutate(genotype=factor(genotype,levels=c('DD','GG','GD')),condition=factor(condition,levels=c('0 h','24 h')))
sm<-iso%>%group_by(genotype,condition)%>%summarise(mean=mean(isoform_fraction),sem=sd(isoform_fraction)/sqrt(n()),.groups='drop')
q6C<-ggplot(sm,aes(genotype,mean,fill=condition))+geom_col(position=position_dodge(.72),width=.65,colour='grey35',linewidth=.3)+geom_errorbar(aes(ymin=mean-sem,ymax=mean+sem),position=position_dodge(.72),width=.13,linewidth=.3)+geom_point(data=iso,aes(y=isoform_fraction),position=position_jitterdodge(jitter.width=.06,dodge.width=.72,seed=12),colour='grey35',size=1.3)+scale_fill_manual(values=c('0 h'='white','24 h'='grey40'),name=NULL)+labs(title='UBOX5 (gene_3219)',x='Genotype (2-year)',y='Isoform fraction')+theme(legend.position='bottom')
change<-sm%>%select(genotype,condition,mean)%>%pivot_wider(names_from=condition,values_from=mean)%>%mutate(difference=`24 h`-`0 h`);mp<-mean(change$difference[change$genotype!='GD'])
q6D<-ggplot(change,aes(genotype,difference,fill=genotype))+geom_col(width=.65,colour='grey40',linewidth=.3)+geom_hline(yintercept=0,colour='grey50')+geom_hline(yintercept=mp,linetype='dashed',colour='grey40')+annotate('text',x=2,y=.38,label=sprintf('Mid-parent ΔIF = %.3f',mp),size=2.8)+annotate('text',x=2,y=-.65,label='ΔΔIF = −0.773',colour='#E41A1C',size=3,fontface='bold')+scale_fill_manual(values=pal,guide='none')+coord_cartesian(ylim=c(-.7,.45))+labs(x='Genotype (2-year)',y='Mean ΔIF (24 h − 0 h)')
apa<-tab('APA_descriptive_values');names(apa)<-make.names(names(apa));print(names(apa))
# 列名从实际表选择，不猜测值向量的顺序。
ac<-names(apa)[grepl('age',names(apa),ignore.case=TRUE)][1];vc<-names(apa)[grepl('delta.*dui|dddui|ΔΔ',names(apa),ignore.case=TRUE)][1]
if(is.na(vc))stop('APA value column not found')
apa$Age<-factor(gsub(' year','-year',as.character(apa[[ac]])),levels=ages)
if(anyNA(apa$Age))apa$Age<-factor(apa[[ac]],levels=c(1,2),labels=ages)
q6E<-ggplot(apa,aes(Age,.data[[vc]],fill=Age))+geom_violin(linewidth=.4,scale='width')+geom_boxplot(width=.17,outlier.shape=NA,fill='white',linewidth=.4)+geom_hline(yintercept=0,linetype='dashed',colour='#E41A1C')+scale_fill_manual(values=c('1-year'='#FDB462','2-year'='#B3B3B3'),guide='none')+labs(x='Age class',y='ΔΔDUI')
ad<-apa%>%mutate(direction=case_when(.data[[vc]]<0~'Negative hybrid deviation',.data[[vc]]>0~'Positive hybrid deviation',TRUE~'No hybrid deviation'))%>%count(Age,direction)%>%group_by(Age)%>%mutate(pct=100*n/sum(n))%>%ungroup()%>%mutate(direction=factor(direction,levels=c('Positive hybrid deviation','No hybrid deviation','Negative hybrid deviation')))
q6F<-ggplot(ad,aes(Age,pct,fill=direction))+geom_col(width=.65,colour='grey40',linewidth=.2)+geom_text(aes(label=ifelse(pct>1,sprintf('%.1f%%',pct),'')),position=position_stack(vjust=.5),size=3)+scale_fill_manual(values=c('Positive hybrid deviation'='#E69F00','Negative hybrid deviation'='#CC79A7','No hybrid deviation'='#B3B3B3'),name=NULL)+labs(x='Age class',y='Proportion of genes (%)')+theme(legend.position='bottom',legend.key.height=unit(2,'mm'))+guides(fill=guide_legend(ncol=1))
fig6<-plot_grid(q6A,q6B,q6C,q6D,q6E,q6F,ncol=2,labels=LETTERS[1:6],label_size=15,label_fontface='plain',hjust=0,vjust=1)
savefig(fig6,6,180,205)

enrich2<-tab('OD_module_enrichment_2yr');en<-bind_rows(mutate(tab('OD_module_enrichment_1yr'),Age='1-year'),mutate(enrich2,Age='2-year'))%>%filter(module%in%c('skyblue','turquoise','darkorange','tan','red','paleturquoise','navajowhite','darkturquoise','white','yellow','blue','green','plum'))%>%mutate(module=factor(module,levels=rev(c('skyblue','turquoise','darkorange','tan','red','paleturquoise','navajowhite','darkturquoise','white','yellow','blue','green','plum'))))
q7A<-ggplot(en,aes(Age,module,size=od_overlap,fill=fold_enrichment,colour=p_adj<.05))+geom_point(shape=21,stroke=.5)+scale_size_area(max_size=8,name='OD-high overlap')+scale_fill_gradient(low='#FDD49E',high='#B30000',name='Fold enrichment')+scale_colour_manual(values=c('TRUE'='black','FALSE'='grey70'),labels=c('No','Yes'),name='FDR < 0.05')+labs(x='Age class',y='WGCNA module')+theme(legend.key.height=unit(3,'mm'))
gm<-tab('WGCNA_1e9_group_means')%>%mutate(genotype=factor(genotype,levels=c('DD','GG','GD')),age=factor(age,levels=1:2,labels=ages),lab=paste(genotype,age))
q7B<-ggplot(gm,aes(plum,CZTB,colour=genotype,shape=age))+geom_point(size=2.8)+geom_text(aes(label=lab),nudge_y=38,size=2.5,show.legend=FALSE)+scale_colour_manual(values=pal,name='Genotype')+scale_shape_manual(values=c(16,17),name='Age')+annotate('label',x=.18,y=1320,label='rho = −0.886\nReference P = 0.033\nAge-blocked P = 0.083\nFDR = 1.00',size=2.5,hjust=1,fill='white',colour='black')+coord_cartesian(ylim=c(620,1400))+labs(x='Plum eigengene (24 h group mean)',y='CZTB (min)')+theme(legend.position='bottom')
q7B<-q7B+scale_x_continuous(limits=c(-.26,.29))+guides(colour=guide_legend(nrow=1),shape=guide_legend(nrow=1))+theme(legend.box='vertical',legend.spacing.y=unit(0,'mm'))
rep<-tab('Representative_expression');q7C<-crop_panel(7,0,.503,1,.497)
fig7<-plot_grid(plot_grid(tag(q7A,'A'),tag(q7B,'B'),ncol=2),q7C,ncol=1,rel_heights=c(1.2,1))
savefig(fig7,7,180,182)

# Figure 8 is generated separately by redraw_qpcr_profiles.R.
cat('Molecular figures 4–7 exported\n');print(warnings())
