"""从随附36个RNA池特征值复算补充表11；不涉及原始测序。"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import f
from data_io import read_sheet

D=Path(__file__).resolve().parents[1]
x=read_sheet(D/'Source_Data_Revised.xlsx',sheet_name='WGCNA_library_eigengenes')
ref=read_sheet(D/'Supplementary_Tables_Revised.xlsx',sheet_name='Table 11',skiprows=3).iloc[:2]
rows=[]
for age,g in x.groupby('age',sort=True):
    g=g.copy();y=g['plum'].to_numpy(float); n=len(g)
    geno=pd.get_dummies(g['genotype'],drop_first=True,dtype=float).to_numpy()
    cond=(g['condition']=='Infected').to_numpy(float).reshape(-1,1)
    reduced=np.column_stack([np.ones(n),geno,cond])
    full=np.column_stack([reduced,geno*cond])
    rss=lambda z: float(np.sum((y-z@np.linalg.lstsq(z,y,rcond=None)[0])**2))
    ra,rf=rss(reduced),rss(full);dfn=np.linalg.matrix_rank(full)-np.linalg.matrix_rank(reduced);dfd=n-np.linalg.matrix_rank(full)
    F=((ra-rf)/dfn)/(rf/dfd);p=float(f.sf(F,dfn,dfd))
    r=ref.loc[ref['age']==age].iloc[0]
    for key,value in [('F',F),('p_nominal',p),('SSE_additive',ra),('SSE_full',rf)]:
        assert np.isclose(float(r[key]),value,rtol=1e-10,atol=1e-14),(age,key,r[key],value)
    rows.append(dict(age=age,n_libraries=n,df_numerator=int(dfn),df_residual=int(dfd),F=F,p_nominal=p,SSE_additive=ra,SSE_full=rf))
out=Path(__file__).resolve().parent/'recomputed';out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out/'plum_interaction.csv',index=False)
(out/'plum_interaction_verification.json').write_text(json.dumps({'values_equal':True,'rows':rows},indent=2),encoding='utf-8')
print(json.dumps(rows,indent=2))
