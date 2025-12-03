import pandas as pd

# Path al file prodotto da LAGKINEMATIC 1.0
path = r"OUT_4DMED/trajectories/lagkin1.0.dat"

# Colonne del file 1.0
cols_1 = [
    "time", "id",
    "X1", "Y1", "Z1",
    "X2", "Y2", "Z2",
    "age"
]

# Legge il file .dat ignorando spazi multipli
df_old = pd.read_csv(
    path,
    sep=r"\s+",
    header=None,
    names=cols_1,
    engine="python"
)

# Il 1.0 mette le depth NEGATIVE: li rendiamo POSITIVE (come nel 3.0)
df_old["Z1"] = -df_old["Z1"]
df_old["Z2"] = -df_old["Z2"]

# Ordina le colonne come nella versione 3.0
cols_3 = ["id", "time", "X1", "Y1", "Z1", "X2", "Y2", "Z2", "age"]
df_old = df_old[cols_3]

# Stampa
print(df_old)

# Salva come CSV
df_old.to_csv("lagkin1_converted.csv", index=False)
print("File salvato: lagkin1_converted.csv")
