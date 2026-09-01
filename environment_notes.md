# Environment notes

These notes summarize runtime information recovered from project records and retained outputs. Package versions are listed only when they were recoverable from existing project records or verified container metadata.

## Reproducibility model

This repository supports layered manuscript-level reproducibility. Lightweight reproduction uses processed summary tables and figure/table scripts. Full upstream reruns require external raw data, cisTarget resources, motif annotation files, Docker images, and additional disk/memory.

## pySCENIC/AUCell runtime

- Docker image: `aertslab/pyscenic:0.12.1`.
- GRN inference: GRNBoost2.
- Workers: 4 for GRN, ctx, and AUCell in the retained discovery run.
- Input expression matrix: 2,322 astrocytes by 36,601 genes.
- AUCell output: matched to 2,322 astrocytes.
- Final downstream regulons: 105 pySCENIC regulons.
- Fixed seed: no fixed `--seed` argument was recovered for the retained GRN or AUCell commands; bitwise-identical reruns should not be assumed.

## pySCENIC resources

The retained run used an AertsLab cisTarget hg38 gene-based ranking database with 10 kb upstream/downstream of transcription start sites and v10 motif annotation.

Full retained resource filenames are listed in `configs/pyscenic_runtime_and_resources.md`. These large third-party resources are not redistributed in this repository.

## CellOracle runtime

- Docker image: `kenjikamimoto126/celloracle_ubuntu:0.18.0`.
- CellOracle version: 0.18.0.
- Verified Python version in the retained CellOracle container: Python 3.10.11.
- Verified container package versions include `celloracle==0.18.0`, `scanpy==1.10.0`, `anndata==0.10.6`, `numpy==1.26.4`, `pandas==1.5.3`, `scipy==1.12.0`, `h5py==3.10.0`, `scikit-learn==1.3.0`, `matplotlib==3.6.3`, `seaborn==0.13.2`, and `statsmodels==0.14.0`.
- CellOracle outputs are in silico perturbation simulations and should be interpreted as candidate-prioritization evidence only.

Approx. Recovery Index / Approx. RI is a study-defined, non-negative composite score calculated as `mean_net_shift_overall * (1 + group_selectivity) * (1 + coherence_target_group) * direction_agreement_factor`, where the direction-agreement factor is 1.0 for concordant group direction and 0.75 otherwise. It is not a cosine similarity, is not normalized to 0-1, and has the theoretical range `[0, +inf)`. It is intended only for within-pipeline candidate comparison and is not a CellOracle official standard metric or a clinical recovery metric.

Formal permutation controls use 10,000 permutations with fixed seed `20240817`. They assess observed-versus-randomized state-shift magnitude and group-specific vector directionality; they do not replace independent biological replication.

## Python package notes

The workflow uses `pandas`, `numpy`, `scipy`, `scanpy`, `anndata`, `h5py`, `matplotlib`, `seaborn`, `scikit-learn`, `celloracle`, `pyscenic`, `statsmodels`, and `openpyxl`.

Version status:

- Recovered from the CellOracle container: `pandas`, `numpy`, `scipy`, `scanpy`, `anndata`, `h5py`, `matplotlib`, `seaborn`, `scikit-learn`, `celloracle`, and `statsmodels`.
- Recovered from the pySCENIC Docker image tag: `pyscenic==0.12.1`.
- `openpyxl`: not pinned / not recovered for the original analysis runtime.

`requirements.lock.txt` records the recovered package pins and marks unrecovered packages explicitly.

## Enrichr / DSigDB / SwissADME

Functional enrichment used Enrichr web API calls for GO/KEGG analyses. Exploratory drug-signature annotation used Enrichr DSigDB results and SwissADME blood-brain barrier annotations. No official Enrichr API version or DSigDB database release number was recovered from retained outputs.

## Excluded resources

The repository intentionally excludes large raw single-cell matrices, controlled-access raw matrices, cisTarget ranking databases, motif annotation files, Docker image layers, and large intermediate objects. See `MANIFEST.md` and `results/audit/repository_release_audit.md` for the release audit.
