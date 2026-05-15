import sys
import csv
import gzip
import numpy as np
from scipy.io import mmread

if len(sys.argv) != 5:
    print("Usage: python mtx_to_csv.py <mtx> <barcodes.tsv> <genes.tsv> <out.csv>", file=sys.stderr)
    sys.exit(1)

mtx_path, barcodes_path, genes_path, out_csv = sys.argv[1:]

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")

def open_binary(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")

def read_barcodes(path):
    vals = []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            vals.append(line.split("\t")[0])
    return vals

def read_genes(path):
    rows = []
    with open_text(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))

    if not rows:
        return []

    ncols = max(len(r) for r in rows)

    # 如果是 2/3 列格式，优先取第 2 列作为 gene symbol
    if ncols >= 2:
        second_col_nonempty = sum(1 for r in rows if len(r) > 1 and r[1].strip())
        if second_col_nonempty >= max(1, len(rows) // 2):
            genes = []
            for r in rows:
                if len(r) > 1 and r[1].strip():
                    genes.append(r[1].strip())
                else:
                    genes.append(r[0].strip() if len(r) > 0 else "NA")
            return genes

    return [r[0].strip() if len(r) > 0 else "NA" for r in rows]

def make_unique(names):
    seen = {}
    out = []
    dup_n = 0
    for x in names:
        x = x.strip() if x else "NA"
        if x not in seen:
            seen[x] = 0
            out.append(x)
        else:
            seen[x] += 1
            dup_n += 1
            out.append(f"{x}__dup{seen[x]}")
    return out, dup_n

print("[1/5] Reading matrix...", file=sys.stderr)
with open_binary(mtx_path) as f:
    mat = mmread(f)

if not hasattr(mat, "tocsr"):
    raise ValueError("Matrix is not sparse or unreadable.")

print("[2/5] Reading barcodes and genes...", file=sys.stderr)
barcodes = read_barcodes(barcodes_path)
genes_raw = read_genes(genes_path)
genes, dup_n = make_unique(genes_raw)

n_rows, n_cols = mat.shape
print(f"Matrix shape: {n_rows} x {n_cols}", file=sys.stderr)
print(f"Barcodes: {len(barcodes)}", file=sys.stderr)
print(f"Genes: {len(genes)}", file=sys.stderr)
print(f"Duplicated gene names fixed: {dup_n}", file=sys.stderr)

# 自动判断方向
# 情况A：rows = cells, cols = genes
if n_rows == len(barcodes) and n_cols == len(genes):
    orientation = "cell_by_gene"
# 情况B：rows = genes, cols = cells
elif n_rows == len(genes) and n_cols == len(barcodes):
    orientation = "gene_by_cell"
else:
    print("ERROR: matrix dimensions do not match barcodes/genes counts.", file=sys.stderr)
    print(f"matrix rows={n_rows}, cols={n_cols}", file=sys.stderr)
    print(f"barcodes={len(barcodes)}, genes={len(genes)}", file=sys.stderr)
    sys.exit(2)

print(f"[3/5] Detected orientation: {orientation}", file=sys.stderr)

print("[4/5] Writing CSV...", file=sys.stderr)
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["CellID"] + genes)

    if orientation == "cell_by_gene":
        mat = mat.tocsr()
        n_cells = n_rows
        n_genes = n_cols
        for i in range(n_cells):
            row = mat.getrow(i)
            dense = np.zeros(n_genes, dtype=np.int64)
            if row.nnz > 0:
                dense[row.indices] = row.data.astype(np.int64, copy=False)
            writer.writerow([barcodes[i]] + dense.tolist())
            if (i + 1) % 100 == 0 or (i + 1) == n_cells:
                print(f"  wrote {i+1}/{n_cells} cells", file=sys.stderr)

    else:
        # gene_by_cell: 按列遍历，每一列是一个 cell
        mat = mat.tocsc()
        n_genes = n_rows
        n_cells = n_cols
        for j in range(n_cells):
            col = mat.getcol(j)
            dense = np.zeros(n_genes, dtype=np.int64)
            if col.nnz > 0:
                dense[col.indices] = col.data.astype(np.int64, copy=False)
            writer.writerow([barcodes[j]] + dense.tolist())
            if (j + 1) % 100 == 0 or (j + 1) == n_cells:
                print(f"  wrote {j+1}/{n_cells} cells", file=sys.stderr)

print("[5/5] Done.", file=sys.stderr)
print(f"Output: {out_csv}", file=sys.stderr)
