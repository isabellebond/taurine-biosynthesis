import pandas as pd

df = pd.read_csv('/home/rmgpibo/Scratch/tmp/n_acetyltaurine_36357675.preds', sep = '\t')

df = df.loc[df['PoPS_Score']> 1]
df = df.sort_values(by = 'PoPS_Score', ascending=False)
print(df.head)