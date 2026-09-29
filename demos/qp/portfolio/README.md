# QP Portfolio Optimization Demo

## Problem

Continuous quadratic optimization applied to portfolio optimization.

## Dataset

MINLPLib `qp3`.

Expected dataset:

`demos/datasets/portfolio/raw/qp3.mps`

## Solver

Repository QP solver.

## Workflow

1. Load the MPS instance.
2. Validate the mathematical model.
3. Verify that the model is continuous.
4. Verify that the objective contains quadratic terms.
5. Translate the model into the repository QP representation.
6. Solve using the QP solver.
7. Report the solution status and solver statistics.

The same portfolio instance is also used by the SuperADMM demo to provide a direct comparison between the two QP solution approaches.
