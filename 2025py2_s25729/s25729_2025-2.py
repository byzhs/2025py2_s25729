from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

def main():
    Entrez.email = input("Email: ")
    Entrez.api_key = input("API key: ")
    t = input("TaxID: ")
    a = int(input("Min length: "))
    b = int(input("Max length: "))
    try:
        h = Entrez.esearch(db="nucleotide", term=f"txid{t}[Organism]", usehistory="y")
        r = Entrez.read(h)
        c = int(r["Count"])
        w, q = r["WebEnv"], r["QueryKey"]
    except Exception as e:
        print("error:", e)
        return
    d, m, L = [], 100, 1000
    for s in range(0, min(c, L), m):
        try:
            h = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=s, retmax=m, webenv=w, query_key=q)
            for x in SeqIO.parse(h, "gb"):
                l = len(x.seq)
                if a <= l <= b:
                    d.append({"Accession": x.id, "Length": l, "Description": x.description})
            h.close()
        except:
            continue
    if not d:
        print("No records"); return
    df = pd.DataFrame(d)
    df.to_csv("records.csv", index=False)
    df = df.sort_values("Length", ascending=False)
    plt.figure(figsize=(12,6))
    plt.plot(df["Accession"], df["Length"], marker='o')
    plt.xticks(range(0, len(df), 10), df["Accession"].iloc[::10], rotation=45)
    plt.tight_layout()
    plt.savefig("sequence.png")

if __name__ == "__main__":
    main()