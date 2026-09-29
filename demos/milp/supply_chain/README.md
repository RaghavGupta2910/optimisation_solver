# MILP Supply Chain Demo

## Problem

Supply-chain and logistics optimization formulated as a mixed-integer linear program.

## Dataset

MIPLIB `shs1023`.

Expected dataset:

`demos/datasets/supply_chain/raw/shs1023.mps`

## Solver

Repository MILP / Branch-and-Cut solver.

## Workflow

1. Load the MPS instance.
2. Validate the mathematical model.
3. Verify that integer variables are present.
4. Solve using the MILP engine.
5. Report the solution status and solver statistics.

The demo does not relax the MILP into an LP.
