# pySCENIC runtime and resource notes

## Discovery run

- Docker image: `aertslab/pyscenic:0.12.1`.
- GRN inference: GRNBoost2.
- Workers: `--num_workers 4` for GRN, ctx, and AUCell.
- Input expression matrix: `input/expression_for_pyscenic.csv`, representing 2,322 astrocytes by 36,601 genes. This large derived matrix is not redistributed.
- AUCell matching: AUCell output was matched back to 2,322 astrocytes with no missing or extra cells in the retained merge record.
- Final downstream regulons: 105.
- Fixed seed: no fixed seed was recovered for the retained GRN or AUCell commands.

## cisTarget / motif resources retained in the run

Repository documentation should describe the retained resources as:

> cisTarget hg38 gene-based ranking database with 10 kb upstream/downstream of transcription start sites and v10 motif annotation.

Full retained filenames:

- TF list: `db/allTFs_hg38.txt`.
- Ranking database: `db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`.
- Motif annotation: `db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`.

These resources are third-party files and are not redistributed in the manuscript-level release package.

## Manuscript consistency note

`Manuscript.docx` has been updated to match the retained v10 resource description. No unresolved pySCENIC resource inconsistency remains in the audited manuscript package.
