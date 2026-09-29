# PDLP Transportation Demo

## Problem

Transportation optimization formulated as a linear program.

## Dataset

Public transportation LP benchmark.

Expected dataset:

`demos/datasets/transportation/raw/transportation.mps`

> The exact public LP instance is still to be verified before adding the dataset.

## Solver

Repository PDLP solver.

## Workflow

1. Load the MPS instance.
2. Validate the mathematical model.
3. Verify that the model contains no integer variables.
4. Verify that the model contains no quadratic objective terms.
5. Solve using the PDLP engine.
6. Report the solution status and solver statistics.

The demo is restricted to a genuine linear-programming formulation and does not silently relax or transform an integer or quadratic model.
