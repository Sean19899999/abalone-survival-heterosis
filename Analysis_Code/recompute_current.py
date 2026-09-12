"""从随附的当前处理数据复算统计；不重建原动物记录或原始组学流水线。"""
from pathlib import Path
import json,itertools
import numpy as np
import pandas as pd
from scipy import stats
from data_io import read_sheet
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent/'recomputed';OUT.mkdir(exist_ok=True)
WB=ROOT/'Source_Data_Revised.xlsx'
tab=lambda name:read_sheet(WB,sheet_name=name)
def check(name,a,b):
 assert list(a.columns)==list(b.columns),(name,'columns')
 assert a.shape==b.shape,(name,a.shape,b.shape)
 for c in a:
  if pd.api.types.is_numeric_dtype(a[c]):
   assert np.allclose(a[c].to_numpy(float),pd.to_numeric(b[c]).to_numpy(float),equal_nan=True,rtol=1e-10,atol=1e-11),(name,c)
  else:assert a[c].fillna('').astype(str).tolist()==b[c].fillna('').astype(str).tolist(),(name,c)
 a.to_csv(OUT/(name+'.csv'),index=False)
 return {'name':name,'rows':len(a),'columns':len(a.columns),'values_equal':True}
report=[]
surv=tab('Current_Survival_Tanks')
surv['AUC_proportion_h']=np.trapezoid(surv[['s0','s24','s48','s72','s96']].to_numpy()/100,x=[0,24,48,72,96],axis=1)
report.append(check('survival_tanks',surv,tab('Current_Survival_Tanks')))
rows=[]
for age,f in surv.groupby('age'):
 groups={g:v.s96.to_numpy() for g,v in f.groupby('genotype')};w={'GD':1,'DD':-.5,'GG':-.5}
 diff=sum(w[g]*np.mean(v) for g,v in groups.items());terms={g:w[g]**2*np.var(v,ddof=1)/len(v) for g,v in groups.items()}
 var=sum(terms.values());df=var**2/sum(terms[g]**2/(len(v)-1) for g,v in groups.items());se=np.sqrt(var);ci=stats.t.ppf(.975,df)*se
 rows.append(dict(age=age,estimate_pp=diff,SE=se,df=df,P_nominal=2*stats.t.sf(abs(diff/se),df),CI95_low=diff-ci,CI95_high=diff+ci,n_tanks_per_genotype=3))
contr=pd.DataFrame(rows);ix=np.argsort(contr.P_nominal);contr['P_Holm_two']=0.;contr.loc[ix,'P_Holm_two']=np.minimum(1,np.maximum.accumulate(contr.P_nominal.to_numpy()[ix]*[2,1]))
report.append(check('survival_contrasts',contr,tab('Current_Survival_Contrasts')))
flow=tab('Current_Flow_Source');baseline=flow[flow.Time.eq(0)].groupby(['Genotype','Age']).THC_recorded.mean()
flow['THC_relative_baseline']=[v/baseline.loc[(g,a)] for g,a,v in zip(flow.Genotype,flow.Age,flow.THC_recorded)]
rows=[]
for name,col,unit in [('THC','THC_relative_baseline','fold of genotype-age baseline mean'),('MOR','MOR_ratio_same_animal','percent'),('PHA','PHA_ratio_same_animal','percent'),('ROS','ROS_recorded','fluorescence arbitrary units')]:
 f=flow[['Genotype','Age','Time','Observation_within_time',col]].rename(columns={col:'value'}).copy()
 if name in ('MOR','PHA'):f['value']*=100
 f['measure']=name;f['unit']=unit;rows.append(f)
long=pd.concat(rows,ignore_index=True);report.append(check('flow_values',long,tab('Current_Flow_Values')))
sm=long.groupby(['Genotype','Age','Time','measure','unit'],as_index=False).value.agg(n='count',mean='mean',sd='std');sm['SEM']=sm.sd/np.sqrt(sm.n)
report.append(check('flow_summary',sm,tab('Current_Flow_Summary')))
auc_rows=[]
for (genotype,age,measure,unit),f in sm.groupby(['Genotype','Age','measure','unit']):
 f=f.sort_values('Time')
 auc_rows.append(dict(genotype=genotype,age=age,measure=measure,AUC=float(np.trapezoid(f['mean'],x=f.Time)),unit=unit+' × h'))
pd.DataFrame(auc_rows).to_csv(OUT/'flow_AUC_descriptive.csv',index=False)
qp=tab('Current_RTqPCR_Observations')
assert not any('CZTB' in c for c in qp.columns)
qc=qp.copy()
qc['ct_mean']=qc[['ct_rep1','ct_rep2','ct_rep3']].mean(axis=1)
qc['ct_sd']=qc[['ct_rep1','ct_rep2','ct_rep3']].std(axis=1,ddof=1)
qc['qc_fail']=qc.ct_sd>.5
qc['delta_ct']=qc.ct_mean-qc.ref_geomean
cal=qc[qc.group.eq('GD-1-year')].groupby('gene').delta_ct.mean()
qc['cal_dct']=qc.gene.map(cal)
qc['delta_delta_ct']=qc.delta_ct-qc.cal_dct
qc['relative_expression']=np.exp2(-qc.delta_delta_ct)
report.append(check('RTqPCR_records',qc,qp))
qs=qc.groupby(['gene','genotype','age','group'],as_index=False).agg(n_retained_records=('relative_expression','size'),mean_expr=('relative_expression','mean'),sd_expr=('relative_expression','std'),median_expr=('relative_expression','median'),mean_ct=('ct_mean','mean'),mean_delta_ct=('delta_ct','mean'),sd_delta_ct=('delta_ct','std'))
qs['sem_expr']=qs.sd_expr/np.sqrt(qs.n_retained_records)
qs['sem_delta_ct']=qs.sd_delta_ct/np.sqrt(qs.n_retained_records)
report.append(check('RTqPCR_group_summaries',qs,read_sheet(ROOT/'Supplementary_Tables_Revised.xlsx',sheet_name='Table 9-GroupMean')))
cardiac=tab('Current_Cardiac_Means')
pri=cardiac[['age','genotype','CBTB_min','CZTB_min']].copy()
z=(pri[['CBTB_min','CZTB_min']]-pri[['CBTB_min','CZTB_min']].mean())/pri[['CBTB_min','CZTB_min']].std(ddof=1)
pri['z_CBTB']=z.CBTB_min;pri['z_CZTB']=z.CZTB_min;pri['PRI_H']=(z.CBTB_min+z.CZTB_min)/np.sqrt(2)
report.append(check('PRI_H_six_group_descriptive',pri,tab('Current_PRI_H_Descriptive')))
variance=(1+np.corrcoef(z.T)[0,1])/2
means=tab('WGCNA_1e9_group_means');expected=tab('WGCNA_1e9_all_tests')
perms=np.array(list(itertools.permutations(range(6))));blocked=np.array([a+tuple(i+3 for i in b) for a in itertools.permutations(range(3)) for b in itertools.permutations(range(3))])
rows=[]
for mod in expected.module.unique():
 x=stats.rankdata(means[mod]);x=(x-x.mean())/np.linalg.norm(x-x.mean())
 for trait in ['Survival_AUC','CBTB','CZTB']:
  y=stats.rankdata(means[trait]);y=(y-y.mean())/np.linalg.norm(y-y.mean());rho=float(x@y)
  rows.append(dict(module=mod,trait=trait,n_groups=6,rho=rho,p_exact_720=np.mean(np.abs(y[perms]@x)>=abs(rho)-1e-12),p_age_blocked_36=np.mean(np.abs(y[blocked]@x)>=abs(rho)-1e-12)))
d=pd.DataFrame(rows)
def bh(p):
 p=np.array(p);ix=np.argsort(p);q=np.empty(len(p));q[ix]=np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
d['fdr_exact']=bh(d.p_exact_720);d['fdr_age_blocked']=bh(d.p_age_blocked_36)
report.append(check('WGCNA_1e9_all_tests',d,expected))
assert len(qp)==350 and qp.stored_sample_id.nunique()==46
assert sm.n.min()==4 and sm.n.max()==5
(OUT/'verification.json').write_text(json.dumps({'checks':report,'stored_qpcr_identifiers':46,'PRI_H_PC1_variance_fraction':float(variance),'cardiac':'provided group summaries only; no inferred animal/tank n','qpcr':'Unchanged stored expression records; no time assignment, biological n reconstruction, ANOVA or cardiac association.','scope':'Current processed-data statistics, not a raw sequencing or instrument reanalysis'},indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
