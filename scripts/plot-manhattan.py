
#Plot metabolite GWAS using skyline

import pandas as pd
import numpy as np
import matplotlib as plt
import math
from skyline import (
    grid,
    examples,
    plots,
    coords,
    figure
)
from matplotlib import gridspec
from matplotlib.backends import backend_agg

#get datasets
metabolites = pd.read_csv('data/gwas-norm/metabolite-map.txt', sep = '\t')
print(metabolites)

dfs = []
df_max = []
df_labels = []

for i, row in metabolites.iterrows():
    #read in data
    chunks = []
    columns = ['chr_name', 'start_pos', 'pvalue']
    dtypes = {'chr_name':str, 'start_pos':int, 'pvalue':float}


    for chunk in pd.read_csv(row['Path'], chunksize=10000, compression='gzip', usecols=columns,  sep = '\t', dtype = dtypes):
        chunks.append(chunk)

    df = pd.concat(chunks)
    #Downsample SNPs for plotting
    sig_df = df.loc[np.exp(-df['pvalue']) < 0.1,:]
    non_sig_df = df.loc[np.exp(-df['pvalue']) > 0.1,:]

    #Save a total of 100000 SNPs for plotting
    num_keep = max(100000 - len(sig_df), 0)
    non_sig_df = non_sig_df.sample(n=num_keep, random_state=42)

    df = pd.concat([sig_df, non_sig_df])

    df_max.append(df.pvalue.max())
    dfs.append(df)
    df_labels.append(metabolites['Metabolite'])

ylim = (0, max(df_max) * 1.02)

#Plot a grid on Manhattan plots
spec = gridspec.GridSpec(4, 2)
b38 = coords.HumanGRCh38()

#Grid object to blot to
dg = grid.BaseGrid(spec, canvas=backend_agg.FigureCanvas)

for idx, i in enumerate(dfs, 1):
    print(idx,i)
    x = dg.plot(plots.manhattan, i, b38, ylim = ylim)
    print(x)
    x[0].set_title(metabolites.loc[idx-1,'Metabolite'])

print(type(dg.figure))
gfig = dg.figure
gfig.savefig('results/figs/manhattan.png', dpi = 300)