"""NCD ratio analysis: biotype of the panel genes (plans/ncd_ratio_exploration.md, section 2).

The analysis uses protein-coding genes only. This script looks up the biotype of the 5001 panel genes and writes
data/_utils/gene_biotype.csv, which ncd_ratio_utils.coding_genes reads.

Source: mygene.info (POST v3/query, scope ensembl.gene, in batches; only the public Ensembl gene IDs of the panel
are sent). It returns the Ensembl biotype of the queried gene ID (field ensembl.type_of_gene, e.g. protein_coding,
lncRNA) and the NCBI gene type (type_of_gene). The Ensembl REST service itself was not reachable on 2026-11-18.
IDs that are not found are looked up by gene symbol.
Rule of the analysis: a gene is kept if its Ensembl biotype is protein_coding.

Run once, from this folder, on a machine with internet access:  python gene_biotype.py
The steps of the analysis read the csv and do not need the network.
Output: data/_utils/gene_biotype.csv   gene, ensembl_id, biotype (Ensembl), ncbi_type, found_by (id / symbol /
                                       not found), source, ensembl_release, lookup_date
"""
import datetime
import json
import urllib.parse
import urllib.request

import anndata
import pandas as pd

import ncd_ratio_utils as nr

SERVER = "https://mygene.info/v3"
BATCH = 500
FIELDS = "symbol,type_of_gene,ensembl.gene,ensembl.type_of_gene"


def query(terms, scopes):
    data = urllib.parse.urlencode({"q": ",".join(terms), "scopes": scopes, "fields": FIELDS, "species": "human"}).encode()
    req = urllib.request.Request(SERVER + "/query", data=data, headers={"User-Agent": "python-urllib"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def ensembl_type(hit, ensembl_id=None):
    """Ensembl biotype of a hit: of the queried gene ID if given, otherwise of its first Ensembl entry."""
    ens = hit.get("ensembl")
    if ens is None:
        return None
    ens = ens if isinstance(ens, list) else [ens]
    for e in ens:
        if ensembl_id is None or e.get("gene") == ensembl_id:
            return e.get("type_of_gene")
    return None


def main():
    var = anndata.read_h5ad(f"../../../data/{nr.SAMPLES[0]}/intermediate_data/adata.h5ad", backed="r").var
    table = pd.DataFrame({"gene": var.index.astype(str), "ensembl_id": var["gene_ids"].astype(str).values})
    with urllib.request.urlopen(SERVER + "/metadata", timeout=60) as r:
        release = json.loads(r.read())["src"]["ensembl"]["version"]
    rows = {}
    ids = table["ensembl_id"].tolist()
    for i in range(0, len(ids), BATCH):
        for hit in query(ids[i:i + BATCH], "ensembl.gene"):
            if hit.get("notfound"):
                continue
            b = ensembl_type(hit, hit["query"])
            if b is not None and hit["query"] not in rows:   # first (best scored) hit of an ID
                rows[hit["query"]] = (b, hit.get("type_of_gene"), "id")
        print(f"looked up {min(i + BATCH, len(ids))} of {len(ids)} ids", flush=True)
    by_id = table["ensembl_id"].map(rows)
    table["biotype"] = [r[0] if isinstance(r, tuple) else None for r in by_id]
    table["ncbi_type"] = [r[1] if isinstance(r, tuple) else None for r in by_id]
    table["found_by"] = [r[2] if isinstance(r, tuple) else None for r in by_id]
    missing = table["biotype"].isna()
    if missing.any():
        best = {}
        for hit in query(table.loc[missing, "gene"].tolist(), "symbol"):
            if not hit.get("notfound") and hit["query"] not in best and ensembl_type(hit) is not None:
                best[hit["query"]] = (ensembl_type(hit), hit.get("type_of_gene"))
        for idx in table.index[missing]:
            if table.at[idx, "gene"] in best:
                table.at[idx, "biotype"], table.at[idx, "ncbi_type"] = best[table.at[idx, "gene"]]
                table.at[idx, "found_by"] = "symbol"
    table["found_by"] = table["found_by"].fillna("not found")
    table["source"] = "mygene.info"
    table["ensembl_release"] = release
    table["lookup_date"] = datetime.date.today().isoformat()
    table.to_csv(nr.UTILS_PATH + "gene_biotype.csv", index=False)
    print(f"Ensembl release {release} (via mygene.info)")
    print(table["biotype"].fillna("not found").value_counts().to_string())
    print(table["found_by"].value_counts().to_string())
    print(pd.crosstab(table["biotype"].fillna("not found"), table["ncbi_type"].fillna("none")).to_string())
    print("not protein_coding:", " ".join(f"{g}({b})" for g, b in zip(table["gene"], table["biotype"]) if b != "protein_coding"))


if __name__ == "__main__":
    main()
