#!/usr/bin/env python
# coding: utf-8

# ## STEP1 - Ridurre a SNP indipendenti (lead SNP via LD pruning)

# In[ ]:


import pandas as pd
import numpy as np
import allel


# In[11]:


# 1) prendi solo i 446 SNP, nella matrice genotipica completa
cand_final = pd.read_csv("SNP_finali_Bonferroni_HWE0.05_ctrl_case.csv")
snps = cand_final["SNP"].tolist()

X_sub = X_var[snps].copy()   # X_var = matrice completa (individui x SNP)
print("X_sub shape:", X_sub.shape)

# 2) considera solo i CONTROLLI per calcolare l'LD (standard)
is_ctrl = (y == 0)
G = X_sub.loc[is_ctrl].values.astype("float32")  # (n_ctrl, 446)

n_samples, n_snps = G.shape
print("G shape (controlli):", G.shape)

# 3) Imputazione NaN col genotipo più frequente per SNP
for j in range(n_snps):
    col = G[:, j]
    mask = np.isnan(col)
    if mask.any():
        vals, counts = np.unique(col[~mask], return_counts=True)
        fill = vals[np.argmax(counts)]
        col[mask] = fill
        G[:, j] = col

# 4) LD pruning con locate_unlinked
G_int = G.astype("int8")
gn = G_int.T   # (n_snps, n_samples)

pruned_mask = allel.locate_unlinked(
    gn,
    size=100,   # finestra
    step=50,
    threshold=0.2  # più basso = più aggressivo
)

print("SNP indipendenti (lead):", pruned_mask.sum())

snps_lead = np.array(snps_446)[pruned_mask]

cand_lead = cand_final[cand_final["SNP"].isin(snps_lead)].copy()
cand_lead = cand_lead.sort_values("pval")

print("Lead SNP dopo LD pruning:", cand_lead.shape[0])
cand_lead.to_csv("SNP_lead_LDpruned.csv", index=False)

