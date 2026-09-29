# Demo Datasets

This directory contains the datasets used by the optimization solver application demos.

## Dataset Categories

- `supply_chain/` — MILP supply-chain instance
- `transportation/` — LP transportation instance
- `portfolio/` — QP portfolio optimization instance
- `unit_commitment/` — MIQP unit-commitment instance
- `reactor_optimization/` — NLP reactor optimization instance

Each dataset directory contains:

- `raw/` — original downloaded dataset
- `processed/` — solver-ready or transformed data
- `metadata.yaml` — dataset source and usage information

Dataset files are kept separate from the demo source code so that the same data can be reused by different solver demonstrations.
