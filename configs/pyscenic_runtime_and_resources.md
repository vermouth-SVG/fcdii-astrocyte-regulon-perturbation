# pySCENIC runtime and resource notes

## Discovery run

- Docker image: `aertslab/pyscenic:0.12.1`.
- Command launch environment: Ubuntu-D WSL using Docker Desktop WSL integration.
- Workers: `--num_workers 4` for GRN, ctx, and AUCell.
- GRN method: GRNBoost2 by pySCENIC 0.12.1 default. The recovered GRN command did not explicitly pass `--method`.
- Fixed seed: no fixed `--seed` was supplied in the recovered retained discovery GRN or AUCell commands.

## Input and output files

- Expression matrix: `input/expression_for_pyscenic.csv`, 2,322 cells x 36,601 genes. This large derived matrix is not redistributed.
- TF list: `allTFs_hg38.txt`, 1,892 TFs.
- Ranking database: `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`.
- Motif annotation: `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`.
- GRN output: `output/grn_adj.tsv`. This large adjacency table is not redistributed.
- Regulon output: `output/regulons.csv`. This small regulon table is included under `results/tables/pyscenic/`.
- AUCell output: `output/auc_mtx.csv`. The full per-cell AUC matrix is not redistributed; summary and matched-object audit outputs are included.

## Resource source notes

- TF list source recovered from project audit: `https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt`.
- Ranking database source recovered from project audit: AertsLab cisTarget hg38/refseq_r80/mc_v10_clust gene-based ranking database matching the filename above.
- The direct motif annotation URL was not recovered from the current project files; the local filename and reuse provenance were recovered.
