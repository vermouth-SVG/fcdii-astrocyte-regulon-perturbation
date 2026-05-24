# CellOracle runtime and resource notes

## Runtime

- Docker image: `kenjikamimoto126/celloracle_ubuntu:0.18.0`.
- CellOracle version: 0.18.0.
- Python version: 3.10.11, recovered from the retained CellOracle container.

## Inputs

The CellOracle in silico perturbation module used:

- Matched astrocyte expression object.
- Candidate transcription factors.
- Cell metadata.
- Low-dimensional embedding.
- Inferred regulatory network.

## Outputs

The retained outputs include:

- Perturbation vectors.
- Mean shift.
- Approx. Recovery Index / Approx. RI.
- Random perturbation controls.
- Per-TF and summary perturbation tables and figures.

Some retained paths include legacy names such as `ko_round1` or `round1_ko_*`. These names are retained for path compatibility and provenance. Current manuscript-facing terminology is `in silico perturbation`.

## Approx. RI interpretation

Approx. RI is a study-defined cosine-similarity-based alignment metric:

`Approx. RI = cos(perturbation vector, lesion-to-internal-control reference vector)`.

It is not a CellOracle official standard metric, not an experimental perturbation result, and not a clinical recovery metric. It should be interpreted only as a computational alignment summary within this study.
