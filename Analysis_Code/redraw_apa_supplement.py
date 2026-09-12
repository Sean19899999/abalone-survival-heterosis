"""保留补图4既有样式，仅将方向图例改为相对亲本的偏离。"""
from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from data_io import read_sheet

C=Path(__file__).resolve().parent;D=C.parent;OUT=C/'figure_outputs';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
raw=pd.read_csv(C/'figure_inputs/APA_DUI_per_library.csv');raw=raw[np.isfinite(raw.DUI)]
apa=read_sheet(D/'Source_Data_Revised.xlsx',sheet_name='APA_descriptive_values')
fig=plt.figure(figsize=(6.6,8.6));gs=fig.add_gridspec(3,2,left=.12,right=.97,bottom=.09,top=.94,hspace=.67,wspace=.4)
samples=[];thresholds=[]
for j,age in enumerate(['1-year','2-year']):
    ax=fig.add_subplot(gs[j,:]); rows=raw[raw.age==age]
    ids=sorted(rows['sample'].unique());values=[rows.loc[rows['sample']==s,'DUI'].to_numpy() for s in ids]
    ax.boxplot(values,positions=np.arange(len(ids)),showfliers=False,widths=.55,medianprops={'color':'black'},boxprops={'color':'#567D8C'},whiskerprops={'color':'#567D8C'},capprops={'color':'#567D8C'})
    ax.set(xticks=np.arange(len(ids)),xticklabels=ids,ylabel='DUI',ylim=(-.04,1.04),title=f'{age}: gene-level distribution within each library')
    ax.tick_params(axis='x',rotation=90);ax.text(-.10,1.08,'AB'[j],transform=ax.transAxes,fontweight='bold',fontsize=12)
    for s,v in zip(ids,values):samples.append(dict(age=age,sample=s,n_genes=len(v),median=np.median(v),q25=np.quantile(v,.25),q75=np.quantile(v,.75)))
for j,age in enumerate(['1-year','2-year']):
    ax=fig.add_subplot(gs[2,j]);v=apa.loc[apa.age==age,'delta_delta_DUI'].to_numpy();xs=np.array([0,.05,.1,.15,.2])
    neg=np.array([(v < -t).sum() for t in xs]);pos=np.array([(v > t).sum() for t in xs])
    for label,y,color in [('ΔΔDUI < −threshold',neg,'#CC79A7'),('ΔΔDUI > threshold',pos,'#E69F00')]:ax.plot(xs,100*y/len(v),'o-',ms=3,label=label,color=color)
    for t,n,p in zip(xs,neg,pos):thresholds.append(dict(age=age,threshold=t,n_genes=len(v),n_negative=int(n),n_within=int(len(v)-n-p),n_positive=int(p),negative_percent=100*n/len(v),positive_percent=100*p/len(v)))
    ax.set(xlabel='Absolute ΔΔDUI threshold',ylabel='Assayed genes (%)',ylim=(0,60),title=age,xticks=xs);ax.legend(frameon=False,fontsize=7.5)
    ax.text(-.23,1.08,'CD'[j],transform=ax.transAxes,fontweight='bold',fontsize=12)
for data,name in [(samples,'APA_library_distributions'),(thresholds,'APA_threshold_counts')]:
    result=pd.DataFrame(data); ref=read_sheet(D/'Source_Data_Revised.xlsx',sheet_name=name)
    # 对当前表格逐值核验，不仅检查图片是否生成。
    result=result.sort_values(list(result.columns[:2])).reset_index(drop=True); ref=ref.sort_values(list(result.columns[:2])).reset_index(drop=True)
    assert result.shape==ref.shape,(name,result.shape,ref.shape)
    for col in result:
        if pd.api.types.is_numeric_dtype(result[col]): assert np.allclose(result[col],ref[col],rtol=1e-10,atol=1e-12),(name,col)
        else: assert result[col].tolist()==ref[col].tolist(),(name,col)
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'Supplementary_Figure_4_revised.{ext}',dpi=400,facecolor='white')
plt.close(fig)
print('Supplementary Figure 4 regenerated; all existing library and threshold summaries match.')
