# SuperADMM Portfolio Optimization Demo

## Problem

Continuous quadratic optimization applied to portfolio optimization.

## Dataset

MINLPLib `qp3`.

Expected dataset:

`demos/datasets/portfolio/raw/qp3.mps`

## Solver

Repository SuperADMM QP solver.

## Workflow

1. Load the MPS instance.
2. Validate the mathematical model.
3. Verify that the model is continuous.
4. Verify that the objective contains quadratic terms.
5. Translate the model into the repository QP representation.
6. Solve using SuperADMM.
7. Record convergence information and solver statistics.

The same portfolio instance is used by the standard QP demo, allowing both QP solution approaches to operate on the same mathematical model.
