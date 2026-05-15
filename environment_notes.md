# Environment notes

These notes summarize the runtime information recovered from the project and the local audit environment. Package versions are not pinned unless they were explicitly recoverable from files or command records.

## Detected audit workstation runtimes

- Windows Python detected during package preparation: Python 3.10.4.
- Ubuntu-D WSL Python detected during package preparation: Python 3.12.3.
- Windows R detected during package preparation: R 4.5.2.

These detected versions describe the audit workstation, not necessarily every original analysis environment.

## pySCENIC/AUCell runtime

- Container image used for the retained discovery pySCENIC outputs: `aertslab/pyscenic:0.12.1`.
- Commands were recovered from the Ubuntu-D WSL shell history and used Docker bind mounts under `<PROJECT_ROOT_WSL>`.
- Docker Desktop with Ubuntu-D WSL integration is supported by local provenance logs and configuration, but no retained per-command daemon event was recovered.
- GRN command used `--num_workers 4`.
- ctx command used `--mode custom_multiprocessing` and `--num_workers 4`.
- AUCell command used `--num_workers 4`.
- The recovered retained discovery commands did not include `--seed`; pySCENIC 0.12.1 CLI defaults to `seed=None` for GRN and AUCell when omitted. Exact bitwise reproducibility should therefore not be assumed for those stochastic steps.

## pySCENIC resources

- TF list: `allTFs_hg38.txt`, AertsLab cisTarget TF list, 1,892 entries.
- Ranking database: `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`.
- Motif annotation: `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`.
- Large ranking and motif resources are not redistributed in this GitHub release. Retrieve them from AertsLab/cisTarget resources before full reruns.

## Key Python packages

The scripts use the following packages across the workflow: `anndata`, `celloracle`, `h5py`, `matplotlib`, `numpy`, `pandas`, `pyscenic`, `scanpy`, `scipy`, `seaborn`, and `statsmodels`. See `requirements.txt` for an unpinned package list.

## R notes

R is not required for the included Python-based figure/table regeneration route unless users extend the workflow. R 4.5.2 was detected on the audit workstation, but no R package lockfile was found.

## Large external data and database requirements

Full computational reproduction requires downloading public raw matrices from GEO and large third-party databases from AertsLab/cisTarget. These files are intentionally excluded from this repository and are listed in `excluded_files_manifest.csv` when present locally.

## Known limitations

- This release emphasizes reproducibility of the manuscript-supporting processed outputs and code audit trail, not one-command rerun of all raw-data processing.
- Some original steps were run through Docker/WSL and are represented by recovered command records and redacted logs.
- Steps without fixed seeds or with multiprocessing may not be bitwise reproducible.
- CellOracle perturbation outputs are in silico simulations and should not be interpreted as experimental validation.
