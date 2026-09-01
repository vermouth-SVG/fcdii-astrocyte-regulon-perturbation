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

Approx. RI is a study-defined composite score:

`Approx. RI = mean_net_shift_overall * (1 + group_selectivity) * (1 + coherence_target_group) * direction_agreement_factor`.

The direction-agreement factor is 1.0 when the expected regulon-effect group matches the higher-expression group and 0.75 otherwise. The score is non-negative and unbounded above (`[0, +inf)`); it is not a cosine similarity and is not normalized to 0-1. It is intended only for within-pipeline candidate comparison. It is not a CellOracle official standard metric, an experimental perturbation result, or a clinical recovery metric.

The formal control analysis uses 10,000 paired-label and group-label permutations with fixed seed `20240817`. These tests assess simulated perturbation structure relative to randomized controls, not biological replication across donors.
