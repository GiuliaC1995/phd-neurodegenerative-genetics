#!/usr/bin/env python
# coding: utf-8

# ## Gene-level analysis / prioritization

# In[23]:


lead = pd.read_csv("SNP_lead_LDpruned_with_supermodule.csv")
lead["SNP"] = lead["SNP"].astype(str).str.strip()

lead_annot = pd.read_csv("SNP_lead_LDpruned_annotated_VEP.csv")
lead_annot["SNP"] = lead_annot["SNP"].astype(str).str.strip()

n_lead = lead["SNP"].nunique()
n_lead_annot = lead_annot["SNP"].nunique()

print("SNP lead totali:", n_lead)
print("SNP lead trovati nel VEP:", n_lead_annot)

missing = sorted(set(lead["SNP"]) - set(lead_annot["SNP"]))
print("SNP lead NON trovati nel VEP:", len(missing))
print(missing[:20])


# In[25]:


lead_annot = pd.read_csv("SNP_lead_LDpruned_annotated_VEP.csv")

# Assicuriamoci che le colonne chiave siano pulite
lead_annot["SNP"] = lead_annot["SNP"].astype(str).str.strip()
lead_annot["SYMBOL"] = lead_annot["SYMBOL"].astype(str).str.strip()
lead_annot["super_module"] = lead_annot["super_module"].astype(str).str.strip()

print("Shape lead_annot:", lead_annot.shape)
print(lead_annot[["SNP", "SYMBOL", "super_module"]].head())


# In[27]:


# Tieni solo righe con un gene "vero"
gene_col = "SYMBOL"

lead_annot_gene = lead_annot[
    (lead_annot[gene_col].notna()) &
    (lead_annot[gene_col] != "") &
    (lead_annot[gene_col] != "-")
].copy()

print("Righe con gene valido:", lead_annot_gene.shape[0])
print("SNP unici con gene:", lead_annot_gene["SNP"].nunique())
print("Geni unici:", lead_annot_gene[gene_col].nunique())


# ### Vedere quali geni hanno più SNP e p più basse

# In[29]:


# assumiamo p_col = "pval" (adatta se si chiama diversamente)
p_col = "pval" 
gene_col = "SYMBOL"

gene_summary = (
    lead_annot_gene
    .groupby(gene_col)
    .agg(
        n_snps=("SNP", "count"),
        min_p=(p_col, "min")
    )
    .reset_index()
    .sort_values("min_p")
)

gene_summary.to_csv("gene_summary_leadSNP_VEP.csv", index=False)
gene_summary.head(20)


# In[31]:


gene_modules = (
    lead_annot_gene
    .groupby(gene_col)["super_module"]
    .apply(lambda x: ",".join(sorted(x.unique())))
    .reset_index()
    .rename(columns={"super_module": "modules"})
)

gene_summary_mod = gene_summary.merge(gene_modules, on=gene_col, how="left")
gene_summary_mod.to_csv("gene_summary_leadSNP_VEP_with_modules.csv", index=False)

gene_summary_mod.head(20)


# In[33]:


gs = gene_summary_mod

gs_2plus = gs[gs["n_snps"] >= 2].sort_values("min_p")
gs_2plus.head(20)


# ### Estrarre geni unici

# In[32]:


gs = gene_summary_mod

genes_all = gs["SYMBOL"].dropna().unique()
print(len(genes_all))

# se vuoi i geni con >=2 SNP:
genes_2plus = gs[gs["n_snps"] >= 2]["SYMBOL"].unique()
print(len(genes_2plus))


# In[ ]:


pd.Series(genes_all).to_csv("genes_all_leadSNP.txt", index=False)
pd.Series(genes_2plus).to_csv("genes_2plus_leadSNP.txt", index=False)

