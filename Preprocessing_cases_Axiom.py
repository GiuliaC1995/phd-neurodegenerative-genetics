#!/usr/bin/env python
# coding: utf-8

# ### Caricamento delle librerie

# In[3]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# # **1️⃣ Lettura e caricamento dei dati grezzi**

# In[65]:


geno_cases = pd.read_csv(
    "cases_from_Axiom.txt", #cambia nome file 
    sep="\t",
    skiprows=5,          # ✅ salta le righe dei metadati
    low_memory=False
)

print("✅ File letto correttamente!")
print("Dimensioni:", geno_cases.shape)
geno_cases.head(2)



# In[67]:


geno_cases = geno_cases.drop(columns=["probeset_id"]) #eliminiamo colonna probeset_id presente negli output Axiom
geno_cases.head()


# In[71]:


def compute_alt(ref, a1, a2): #funzione per creare la colonna Alt_Allele se non è presente
    if pd.isna(ref) or pd.isna(a1) or pd.isna(a2):
        return np.nan
    if ref == a1:
        return a2
    elif ref == a2:
        return a1
    else:
        # fallback: assume A2 sia l'alternativo
        return a2

geno_cases["Alt_Allele"] = [
    compute_alt(ref, a1, a2)
    for ref, a1, a2 in zip(geno_cases["Ref_Allele"], geno_cases["Allele_A"], geno_cases["Allele_B"])
]

geno_cases[["Allele_A", "Allele_B", "Ref_Allele", "Alt_Allele"]].head()


# In[73]:


geno_cases.columns.tolist() #controlliamo le colonne del dataset per vedere quale è l'ultima colonna campione


# In[75]:


# per esempio se la prima colonna descrittiva è "Allele_A" 
first_meta_col = "Allele_A"

geno_data = geno_cases.loc[:, :"PD0145.CEL_call_code"] # <-- replace "PD" with your last sample column
annot_data = geno_cases.loc[:, first_meta_col:]

print("Genotipi:", geno_data.shape)
print("Annotazioni:", annot_data.shape)


# In[76]:


annot_reduced = annot_data[["Allele_A", "Allele_B", "Ref_Allele", "Alt_Allele"]]


# In[81]:


geno_cases["BestandRecommended"].value_counts(dropna=False) #vediamo quante varianti hanno BestandRecommended = 1


# In[83]:


geno_cases["ConversionType"].value_counts(dropna=False) #vediamo quante varianti per classe ConversionType


# In[57]:


geno_filtered = geno_cases[
    (geno_cases["BestandRecommended"] == 1) &
    (geno_cases["ConversionType"].isin(["PolyHighResolution"])) #filtriamo tenendo solo varianti PolyHighResolution e BestandRecommended
]

print(f"✅ Varianti mantenute dopo filtro: {geno_filtered.shape[0]}")



# In[19]:


geno_filtered.head()


# In[59]:


# --- Soglie ---
maf_threshold = 0.01
hwe_threshold = 0.05 #cambiare a seconda del disegno di studio

# --- Rimuoviamo varianti con valori mancanti ---
geno_filtered = geno_filtered.dropna(subset=["MinorAlleleFrequency", "H.W.p-Value"])

# --- Applichiamo i filtri MAF e HWE ---
mask_qc = (
    (geno_filtered["MinorAlleleFrequency"] >= maf_threshold) &
    (geno_filtered["H.W.p-Value"] >= hwe_threshold)
)

geno_qc = geno_filtered[mask_qc].copy()

print(f"✅ Varianti dopo filtro MAF ≥ {maf_threshold} e HWE ≥ {hwe_threshold}: {geno_qc.shape[0]}")


# ### Qualche plot

# In[29]:


plt.figure(figsize=(6,4)) #Grafico a barre distribuzione MAF prima del filtro
geno_filtered["MinorAlleleFrequency"].hist(bins=50)
plt.title("Distribuzione MAF (prima del filtro)")
plt.xlabel("Minor Allele Frequency")
plt.ylabel("Conteggio varianti")
plt.show()

plt.figure(figsize=(6,4)) #Grafico a barre distribuzione MAF dopo il filtro
geno_qc["MinorAlleleFrequency"].hist(bins=50)
plt.title("Distribuzione MAF (dopo il filtro)")
plt.xlabel("Minor Allele Frequency")
plt.ylabel("Conteggio varianti")
plt.show()


# In[31]:


plt.figure(figsize=(7,5)) #Boxplot distribuzione MAF per ConversionType 
sns.boxplot(
    data=geno_qc,
    x="ConversionType",
    y="MinorAlleleFrequency",
    palette=["#66c2a5", "#fc8d62"]
)
plt.title("Distribuzione MAF per tipo di ConversionType")
plt.xlabel("ConversionType")
plt.ylabel("Minor Allele Frequency")
plt.show()


# In[23]:


geno_qc["super_module"].value_counts(dropna=False) #vediamo la distribuzione dei vari moduli Axiom


# In[23]:


geno_qc_plot = geno_qc.dropna(subset=["super_module"]).copy() #Plot
geno_qc_plot["super_module"] = geno_qc_plot["super_module"].astype(str)
geno_qc_plot = geno_qc_plot.reset_index(drop=True)


# In[46]:


plt.figure(figsize=(8,5)) #altri plot
sns.countplot(
    data=geno_qc_plot,
    y="super_module",
    order=geno_qc_plot["super_module"].value_counts().index,
    palette="viridis"
)
plt.title("Distribuzione delle varianti per super_module (post QC)")
plt.xlabel("Conteggio varianti")
plt.ylabel("super_module")
plt.show()


# In[48]:


module_counts = (
    geno_qc["super_module"]
    .value_counts(normalize=True)
    .mul(100)
    .round(2)
)

print("📊 Percentuale di varianti per modulo:")
print(module_counts)


# In[25]:


keep_modules = [  #Modificabile a seconda del disegno di studio
    "GWAS Grid",
    "Immunity, Inflammation, and HLA",
    "Published GWAS Hits",
    "Pharmacogenetics/ADME",
    "Wellness and Lifestyle"
]

geno_biomarker = geno_qc[geno_qc["super_module"].isin(keep_modules)].copy()

print(f"✅ Varianti mantenute per analisi biomarcatori: {geno_biomarker.shape[0]:,}")
geno_biomarker["super_module"].value_counts()


# In[52]:


import matplotlib.pyplot as plt

counts = geno_biomarker["super_module"].value_counts(normalize=True).mul(100)
plt.figure(figsize=(7,5))
plt.barh(counts.index, counts.values, color="teal")
plt.xlabel("% varianti")
plt.title("Distribuzione moduli (subset biomarcatori ND)")
plt.show()


# In[29]:


# filtriamo con i supermoduli scelti
geno_stats= geno_qc[geno_qc["super_module"].isin(keep_modules)].copy()

print(f"🧠 Varianti dopo filtro supermodulo: {geno_stats.shape[0]:,} varianti")


# In[31]:


geno_stats.head(10)


# In[33]:


geno_stats = geno_stats.reset_index() #eliminiamo indice
geno_stats.head()


# In[45]:


# === 🔹 Crea dataset allelico mantenendo subito le annotazioni utili ===

# Copia dataset di lavoro
geno = geno_stats.copy()

# Colonne dei genotipi
sample_cols = [c for c in geno.columns if c.startswith("PD")] # <-- replace "PD" with your sample ID prefix

# Funzione per convertire AA/AB/BB → alleli reali
def convert_genotypes(row):
    a = row["Allele_A"]
    b = row["Allele_B"]
    mapping = {
        "AA": a + a,
        "AB": a + b,
        "BB": b + b,
        "NoCall": np.nan,
        "NC": np.nan,
    }
    return [mapping.get(row[col], np.nan) for col in sample_cols]

# Conversione allelica
geno_allelic = geno.apply(
    lambda r: pd.Series(convert_genotypes(r), index=sample_cols),
    axis=1
)

# ===  Mantieniamo ora TUTTE le colonne essenziali ===
geno_allelic.insert(0, "dbSNP_RS_ID", geno["dbSNP_RS_ID"])
geno_allelic.insert(1, "Ref_Allele", geno["Ref_Allele"])
geno_allelic.insert(2, "Alt_Allele", geno["Alt_Allele"])

print("✨ Conversione completa e annotazioni preservate!")
geno_allelic.head()


# # **7️⃣ Controllo e filtro dei valori mancanti**

# In[47]:


# Calcola la percentuale di missing per variante
missing_per_variant = geno_allelic[sample_cols].isna().mean(axis=1) * 100

# Statistiche generali
print("📊 Distribuzione percentuale di missing per SNP:")
print(missing_per_variant.describe())

# Grafico
import matplotlib.pyplot as plt
plt.hist(missing_per_variant, bins=50)
plt.xlabel("% missing per variante")
plt.ylabel("Numero di varianti")
plt.title("Distribuzione missing values per variante")
plt.show()


# In[43]:


# Numero di missing per SNP
missing_counts = geno_allelic[sample_cols].isna().sum(axis=1)

# Statistiche generali
print("📊 Distribuzione numero di missing per SNP:")
print(missing_counts.describe())

# Istogramma (in valori assoluti)
import matplotlib.pyplot as plt
plt.hist(missing_counts, bins=50)
plt.xlabel("Numero di missing per variante")
plt.ylabel("Numero di varianti")
plt.title("Distribuzione missing values per variante (valori assoluti)")
plt.show()


# In[49]:


# === 🔹 Filtro varianti con >5% di missing (NoCall inclusi) ===

threshold = 0.05  # 5% di missing consentito per variante

# Calcola percentuale di missing per ciascuna variante
missing_per_variant = geno_allelic[sample_cols].isna().mean(axis=1)

# Filtra: tieni solo varianti con meno del 5% di missing
geno_allelic_filtered = geno_allelic[missing_per_variant < threshold].copy()

print(f"✅ Varianti dopo rimozione di quelle con >{threshold*100:.0f}% missing: {geno_allelic_filtered.shape[0]}")

# (Facoltativo) Verifica la distribuzione dei missing residui
import matplotlib.pyplot as plt
plt.hist(missing_per_variant[missing_per_variant < threshold] * 100, bins=50)
plt.xlabel("% missing per variante")
plt.ylabel("Numero di varianti")
plt.title("Distribuzione missing values (dopo filtro 5%)")
plt.show()


# In[111]:


# === 🔹 Report finale sui valori mancanti (post-QC) ===

# Percentuale di varianti con almeno un genotipo mancante
variants_with_missing = (
    geno_allelic_filtered[sample_cols].isna().sum(axis=1) > 0
).mean() * 100

# Numero totale di valori mancanti
missing_total = geno_allelic_filtered[sample_cols].isna().sum().sum()

# Percentuale di genotipi mancanti per campione
missing_per_sample = geno_allelic_filtered[sample_cols].isna().sum()
missing_percent = (missing_per_sample / geno_allelic_filtered.shape[0]) * 100

print("📊 --- REPORT QC VALORI MANCANTI ---")
print(f"🔹 Missing totali nel dataset: {missing_total:,}")
print(f"🔹 Varianti con ≥1 missing: {variants_with_missing:.2f}%")
print(f"🔹 Percentuale media di missing per campione: {missing_percent.mean():.2f}%")
print("🔹 Range missing per campione:")
print(f"   Min: {missing_percent.min():.2f}% | Max: {missing_percent.max():.2f}%")

# (Facoltativo) Grafico distribuzione % missing per campione
plt.figure(figsize=(6,4))
plt.hist(missing_percent, bins=40, color="teal", edgecolor="black")
plt.xlabel("% genotipi mancanti per campione")
plt.ylabel("Numero di campioni")
plt.title("Distribuzione della percentuale di genotipi mancanti per campione")
plt.show()


# # **8️⃣ Salvataggio finale dei dataset**

# In[45]:


# === SALVATAGGIO FINALE DEI DATASET QC ===

# 🔹 1️⃣ Dataset completi post-QC (prima della conversione allelica)
geno_stats.to_csv("PD_stats_QC.txt", sep="\t", index=False) # <-- replace with your ouput file

# 🔹 2️⃣ Dataset convertito in alleli (già AA/AB/BB → basi reali)
geno_allelic.to_csv("PD_QC_allelic.txt", sep="\t", index=False) # <-- replace with your ouput file

# 🔹 3️⃣ (Facoltativo) Versioni .pkl per caricamento rapido
geno_stats.to_pickle("PD_stats_QC.pkl") # <-- replace with your ouput file
geno_allelic.to_pickle("PD_stats_QC_allelic.pkl") # <-- replace with your ouput file

print("✅ Tutti i dataset salvati con successo!")
print("📁 File creati:")
print(" - PD_stats_QC_allelic.txt / .pkl") # <-- replace with your ouput file


# ### Conversione genotipi in 0,1,2

# In[3]:


import pandas as pd

geno_allelic_filtered = pd.read_csv(
    "PD_QC_allelic.txt", # <-- replace with your input file
    sep="\t",
)


# In[5]:


geno_allelic_filtered= geno_allelic_filtered.set_index("dbSNP_RS_ID")
geno_allelic_filtered.head()


# In[11]:


geno_allelic_filtered= geno_allelic_filtered.reset_index()
geno_allelic_filtered.head()


# In[43]:


# === 🔹 Ripristina Ref_Allele e Alt_Allele in geno_imputed ===

geno_allelic_filtered = geno_allelic_filtered.merge(
    geno_filtered[["dbSNP_RS_ID", "Ref_Allele", "Alt_Allele"]],
    on="dbSNP_RS_ID",
    how="left"
)

print("✅ Colonne Ref_Allele e Alt_Allele reintegrate!")
print([c for c in geno_allelic_filtered.columns if "Allele" in c])


# In[53]:


# Controlla che non ci siano Ref/Alt mancanti
missing_refalt = geno_allelic_filtered[geno_allelic_filtered["Ref_Allele"].isna()]
print(f"Varianti senza Ref/Alt trovate: {missing_refalt.shape[0]}")


# In[7]:


# === 🔹 Conversione dei genotipi allelici in formato numerico (0/1/2) ===

geno_numeric = geno_allelic_filtered

# Identifica le colonne dei campioni
sample_cols = [c for c in geno_numeric.columns if c.startswith("PD")] #cambiare

# Funzione per codificare ogni riga (variante)
def encode_genotype(row):
    ref = row["Ref_Allele"]
    alt = row["Alt_Allele"]
    mapping = {
        ref + ref: 0,   # omozigote reference
        ref + alt: 1,   # eterozigote
        alt + ref: 1,   # eterozigote inverso
        alt + alt: 2,   # omozigote alternativo
    }
    return [mapping.get(row[col], np.nan) for col in sample_cols]

# Applica la codifica a ogni variante (riga)
geno_encoded = geno_numeric.apply(lambda r: pd.Series(encode_genotype(r), index=sample_cols), axis=1)

# Aggiungi l'ID della variante e le info base
geno_encoded.insert(0, "dbSNP_RS_ID", geno_numeric["dbSNP_RS_ID"])
geno_encoded["Ref_Allele"] = geno_numeric["Ref_Allele"]
geno_encoded["Alt_Allele"] = geno_numeric["Alt_Allele"]

print("✅ Conversione completata!")
print("Esempio:")
print(geno_encoded.head(5))

# --- Controllo rapido ---
print("\n📊 Valori unici nei genotipi numerici:")
print(pd.Series(np.concatenate([geno_encoded[c].dropna().unique() for c in sample_cols])).value_counts())


# In[17]:


geno_encoded.to_pickle("PD_QC_encoded.pkl") # <-- replace with your input file


# In[11]:


# === 🔹 Conta i NaN nei genotipi numerici ===

# Colonne dei campioni
sample_cols = [c for c in geno_encoded.columns if c.startswith("PD")] #cambiare

print("📌 Numero totale di NaN nei genotipi:")
total_nan = geno_encoded[sample_cols].isna().sum().sum()
print(total_nan)

print("\n📌 Percentuale di NaN rispetto a tutto il dataset genotipico:")
total_vals = geno_encoded[sample_cols].size
print((total_nan / total_vals) * 100, "%")

print("\n📌 NaN per variante (prime 5 righe):")
display( geno_encoded[sample_cols].isna().sum(axis=1).head() )

print("\n📌 NaN per campione (prime 5 colonne):")
display( geno_encoded[sample_cols].isna().sum().head() )


