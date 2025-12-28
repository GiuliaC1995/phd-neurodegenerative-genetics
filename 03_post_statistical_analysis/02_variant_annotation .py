#!/usr/bin/env python
# coding: utf-8

# ## STEP 1 - Annotazione Axiom (super_module)

# In[ ]:


annot =pd.read_csv(   # ricarichiamo il file con annotazioni
    "genotyping data.txt",
    sep="\t",
    skiprows=5
)


# In[82]:


lead = pd.read_csv("SNP_lead_LDpruned.csv") # file ottenuto dopo LD pruning
print(lead.shape)
print(lead.columns)
lead.head()


# In[7]:


"super_module" in annot.columns


# In[9]:


# Pulisci i nomi delle colonne
annot.columns = annot.columns.str.strip()
print(annot.columns[:15])

# Tieni solo dbsnp_rs_id e super_module
annot = annot[["dbSNP_RS_ID", "super_module"]].copy()

annot["dbSNP_RS_ID"] = annot["dbSNP_RS_ID"].astype(str).str.strip()
annot["super_module"] = annot["super_module"].astype(str).str.strip()

annot = annot.dropna(subset=["dbSNP_RS_ID", "super_module"])
annot = annot.drop_duplicates(subset=["dbSNP_RS_ID"])

print("Annotazione Axiom:", annot.shape)
annot.head()


# In[11]:


# Dizionario rsID -> super_module
snp_to_module = annot.set_index("dbSNP_RS_ID")["super_module"].to_dict()

# Aggiungi la colonna super_module ai lead
lead["SNP"] = lead["SNP"].astype(str).str.strip()
lead["super_module"] = lead["SNP"].map(snp_to_module)

print("Lead con super_module assegnato:", lead["super_module"].notna().sum())
print("Lead senza modulo (NaN):", lead["super_module"].isna().sum())

# Distribuzione dei moduli
print(lead["super_module"].value_counts(dropna=False))

# Salva
lead.to_csv("SNP_lead_LDpruned_with_supermodule.csv", index=False)


# In[17]:


lead.head()


# ## STEP 2 - ANNOTAZIONE FUNZIONALE

# In[19]:


# 1) Lead SNP già con super_module
lead = pd.read_csv("SNP_lead_LDpruned_with_supermodule.csv")
lead["SNP"] = lead["SNP"].astype(str).str.strip()
print("Lead:", lead.shape)
print(lead.columns)

# 2) VEP
vep = pd.read_csv("vep.txt", sep="\t") # <-- modificare con file ottenuto da VEP
print("Colonne VEP:", vep.columns[:15])

# scegli colonna rsID (Uploaded_variation) in modo robusto
if "#Uploaded_variation" in vep.columns:
    snp_col = "#Uploaded_variation"
elif "Uploaded_variation" in vep.columns:
    snp_col = "Uploaded_variation"
else:
    raise ValueError("Non trovo la colonna Uploaded_variation in VEP")

# scegli alcune colonne utili (adatta se ne hai di diverse)
cols_keep = [snp_col, "SYMBOL"]
for extra in ["Consequence", "BIOTYPE"]:
    if extra in vep.columns:
        cols_keep.append(extra)

vep_sub = (
    vep[cols_keep]
    .dropna(subset=[snp_col, "SYMBOL"])
    .drop_duplicates()
    .rename(columns={snp_col: "SNP"})
)

vep_sub["SNP"] = vep_sub["SNP"].astype(str).str.strip()
vep_sub["SYMBOL"] = vep_sub["SYMBOL"].astype(str).str.strip()

print("VEP sub:", vep_sub.shape)
vep_sub.head()

# 3) Merge: lead + VEP (gene, consequence, biotype)
lead_annot = lead.merge(vep_sub, on="SNP", how="left")

print("Lead annotati con VEP:", lead_annot.shape)
print(lead_annot[["SNP", "super_module", "SYMBOL"]].head())

lead_annot.to_csv("SNP_lead_LDpruned_annotated_VEP.csv", index=False)

