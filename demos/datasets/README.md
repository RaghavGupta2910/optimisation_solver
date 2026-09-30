# Demo Datasets

This directory contains the datasets used by the optimization solver application demos.

## Dataset Categories

- supply_chain/ - MILP supply-chain instance
- transportation/ - LP transportation instance
- portfolio/ - QP portfolio optimization instance
- unit_commitment/ - MIQP unit-commitment instance
- reactor_optimization/ - NLP reactor optimization instance

Each dataset directory contains:

- raw/ - original downloaded dataset
- processed/ - solver-ready or transformed data
- metadata.yaml - dataset source and usage information

Dataset files are kept separate from the demo source code so that the same data can be reused by different solver demonstrations.

## Supply Chain

The supply-chain demo uses the MIPLIB shs1023 benchmark.

The benchmark is included at:

demos/datasets/supply_chain/raw/shs1023.mps

Source:

https://miplib.zib.de/WebData/instances/shs1023.mps.gz

The extracted benchmark is approximately 89 MB.

Run the demo with:

.\build\demos\Debug\supply_chain_milp.exe "demos\datasets\supply_chain\raw\shs1023.mps"

## Transportation

The transportation demo uses the small Dantzig transportation instance:

demos/datasets/transportation/raw/mpstrans.mps

## Portfolio

The portfolio demo currently includes the MINLPLib qp3 input:

demos/datasets/portfolio/raw/qp3.lp

qp3 has an indefinite quadratic objective and is therefore not compatible with the repository's convex-QP/SuperADMM requirements. The QP demo is retained as an integration example but is currently not expected to solve this instance successfully.

## Unsupported Problem Classes

The MIQP and NLP demos document the current solver integration boundaries for those problem classes. They do not silently relax or convert the benchmark into another problem class.
