#!/usr/bin/env python
# coding: utf-8

# ## STEP 1 – Carica casi e controlli, allinea le varianti comuni

# In[3]:


import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from scipy.stats import chisquare


# In[3]:


cases = pd.read_pickle("PD_QC_encoded.pkl") # <-- cambiare con il nome del file ottenuto dal pre-processing
cases.head()


# In[15]:


geno_t = cases.T.copy() # trasponiamo righe e colonne per uniformare al fine dei controlli


# In[17]:


geno_t.head()


# In[19]:


geno_t.columns = geno_t.columns.astype(str) #Converte tutti i nomi delle colonne di geno_t in stringhe (str).


# In[23]:


nan_cols = [ # Remove columns with missing or invalid variant IDs generated after transposition
    c for c in geno_t.columns
    if pd.isna(c) or (isinstance(c, str) and c.lower() == "nan")
]

print(len(nan_cols), nan_cols[:10])

geno_t = geno_t.drop(columns=nan_cols)



# In[24]:


geno_t.shape


# In[27]:


[c for c in geno_t.columns if c.lower().startswith("sample")] #controlliamo se esiste già la colonna 'sample'


# In[29]:


geno_t = geno_t.reset_index() # resettiamo l'indice 

geno_t = geno_t.rename(columns={"index": "Sample"}) # poi rinominiamo la colonna 'index' come 'Sample'

geno_t.head()


# In[35]:


cases=geno_t.copy() # rinominiamo il dataframe


# In[36]:


cases.insert(0, "phenotype", 1) #inseriamo la colonna fenotipo = 1


# In[37]:


cases.head()


# In[41]:


cases = cases.set_index("Sample") # ora impostiamo 'Sample' come indice
cases.head()


# In[42]:


controls = pd.read_pickle("controls_clean.pkl") # <-- cambiare con il nome del file dei controlli
controls.head()


# In[45]:


print("Casi:", cases.shape)
print("Controlli:", controls.shape)

# varianti comuni (escludi phenotype / Sample)
case_vars = set(cases.columns) - {"phenotype", "Sample"}
ctrl_vars = set(controls.columns) - {"phenotype", "Sample"}

common = sorted(list(case_vars & ctrl_vars))
print("Varianti comuni:", len(common))

# allinea e concatena
cases_aligned = cases[["phenotype"] + common]
controls_aligned = controls[["phenotype"] + common]

X = pd.concat([cases_aligned, controls_aligned], ignore_index=True)
print("Nuova X:", X.shape)

# opzionale: salva
X.to_pickle("X_common_variants.pkl") 


# ## STEP 2 – Definisci y e X_var

# In[5]:


X = pd.read_pickle("X_common_variants.pkl")

y = X["phenotype"].values.astype(int)
X_var = X.drop(columns=["phenotype"])
print("X_var shape:", X_var.shape)
print("Distribuzione y:", np.bincount(y))


# ## STEP 3 – TEST DI ASSOCIAZIONE: Fisher genome-wide

# In[13]:


pvals = []

cases_idx = (y == 1)
controls_idx = (y == 0)

for snp in X_var.columns:
    geno = X_var[snp].values  # genotipi 0/1/2

    # normalizza: converte in float → NA → -1 → int
    geno = pd.to_numeric(geno, errors="coerce")
    geno = np.nan_to_num(geno, nan=-1).astype(int)

    # maschera dei validi
    valid_mask = geno >= 0

    geno = geno[valid_mask]
    case_mask = cases_idx[valid_mask]
    ctrl_mask = controls_idx[valid_mask]

    # conteggio allelici
    counts = np.bincount(geno, minlength=3)
    minor = np.argmin(counts)

    case_minor = np.sum((geno == minor) & case_mask)
    case_other = np.sum((geno != minor) & case_mask)

    ctrl_minor = np.sum((geno == minor) & ctrl_mask)
    ctrl_other = np.sum((geno != minor) & ctrl_mask)

    table = np.array([
        [case_minor, case_other],
        [ctrl_minor, ctrl_other]
    ])

    _, p = fisher_exact(table)
    pvals.append(p)

pvals = np.array(pvals)

# FDR Benjamini–Hochberg
_, pvals_adj, _, _ = multipletests(pvals, method='fdr_bh')

results = pd.DataFrame({
    "SNP": X_var.columns,
    "pval": pvals,
    "padj": pvals_adj
}).sort_values("padj")

print(results.head())

# salva risultati completi
results.to_excel("fisher_all_snps.xlsx", index=False)

# SNP significativi a FDR 5%
sig_snps = results[results["padj"] < 0.05]["SNP"].tolist()
print("Varianti significative FDR<0.05:", len(sig_snps))
    


# In[14]:


X_var.to_pickle("X_var")
results[results["padj"] < 0.05].to_excel("fisher_significant.xlsx", index=False)


# ## STEP 4 – QC addizionale (Bonferroni + HWE) sui candidati di Fisher

# ### 4.1 – Carica i candidati da Fisher e applica Bonferroni

# In[19]:


# Carico le varianti significative a FDR<0.05 con Fisher
fisher_sig = pd.read_excel("fisher_significant.xlsx")
print("fisher_sig shape:", fisher_sig.shape)
print(fisher_sig.head())

# fisher_sig contiene almeno: "SNP", "pval", "padj"
# Per lo STEP 5 usiamo il p-value grezzo di Fisher come "pval"
cand_df = fisher_sig.rename(columns={"padj": "p_fdr_fisher"}).copy()

m = cand_df.shape[0]   # numero di SNP candidati
print("Numero di SNP candidati (Fisher FDR<0.05):", m)

# Bonferroni sui pval di Fisher
cand_df["p_bonf"] = np.minimum(cand_df["pval"] * m, 1.0)
cand_df["bonf_signif"] = cand_df["p_bonf"] < 0.05

print("SNP significativi dopo Bonferroni (p_bonf < 0.05):",
      cand_df["bonf_signif"].sum())  # Bonferroni correction applied only to Fisher-significant SNPs (post hoc stringency)


# ### 4.2 – Funzione per HWE (nei controlli e nei casi)

# In[21]:


def hwe_chisq_pval(genotypes):
    """
    Calcola il p-value HWE (chi-quadro 1 df) per genotipi codificati 0/1/2.
    Rimuove i NaN, ritorna np.nan se non calcolabile.
    """
    g = genotypes[~np.isnan(genotypes)]
    if g.size == 0:
        return np.nan

    nAA = np.sum(g == 0)
    nAB = np.sum(g == 1)
    nBB = np.sum(g == 2)
    n = nAA + nAB + nBB
    if n == 0:
        return np.nan

    # frequenza allelica (allele "A" = 0, "B" = 2)
    p = (2 * nAA + nAB) / (2.0 * n)
    q = 1.0 - p

    # attesi sotto HWE
    expAA = n * p * p
    expAB = 2 * n * p * q
    expBB = n * q * q

    obs = np.array([nAA, nAB, nBB], dtype=float)
    exp = np.array([expAA, expAB, expBB], dtype=float)

    if np.any(exp == 0):
        return np.nan

    chi2, pval = chisquare(f_obs=obs, f_exp=exp, ddof=1)
    return pval


# ### 4.3 – Calcola HWE p-value nei controlli e nei casi

# In[23]:


is_ctrl = (y == 0)
is_case = (y == 1)

hwe_ctrl_p = []
hwe_case_p = []
missing_snps = []

for snp in cand_df["SNP"]:
    if snp not in X_var.columns:
        hwe_ctrl_p.append(np.nan)
        hwe_case_p.append(np.nan)
        missing_snps.append(snp)
        continue

    g = X_var[snp].values.astype(float)

    hwe_ctrl_p.append(hwe_chisq_pval(g[is_ctrl]))
    hwe_case_p.append(hwe_chisq_pval(g[is_case]))

cand_df["HWE_p_ctrl"] = hwe_ctrl_p
cand_df["HWE_p_case"] = hwe_case_p

print("SNP non trovati in X_var:", len(missing_snps))


# ### 4.4 – Applica filtri HWE nei controlli e nei casi + Bonferroni

# In[ ]:


alpha_hwe = 0.05

mask_qc_b = (
    cand_df["bonf_signif"] &
    (cand_df["HWE_p_ctrl"] >= alpha_hwe) &
    (cand_df["HWE_p_case"] >= alpha_hwe)
)

cand_final_b = cand_df[mask_qc_b].copy()

print("SNP finali (Bonferroni + HWE_ctrl+case>=0.05):",
      cand_final_b.shape[0])

cand_final_b.to_csv("SNP_finali_Bonferroni_HWE0.05_ctrl_case.csv", index=False)

