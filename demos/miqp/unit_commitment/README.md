# MIQP Unit Commitment Demo

## Problem

Unit-commitment optimization formulated as a mixed-integer quadratic program.

## Dataset

Public MIQP unit-commitment benchmark.

The exact public instance is still to be verified before adding the dataset.

Expected dataset location:

`demos/datasets/unit_commitment/raw/`

## Solver

Repository unified solver interface.

## Workflow

1. Load the MIQP model when a verified public instance is available.
2. Validate the mathematical model.
3. Verify that integer variables are present.
4. Verify that quadratic objective terms are present.
5. Pass the model to the unified solver interface.
6. Report the resulting solver status.

The demo does not silently relax the MIQP into an LP or continuous QP.

If the current solver stack does not provide a complete MIQP branch-and-cut backend, the demo reports that backend boundary explicitly rather than presenting a relaxed problem as an MIQP solution.
