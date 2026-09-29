# NLP Reactor Optimization Demo

## Problem

Nonlinear reactor/process optimization benchmark.

## Dataset

MINLPLib `ex8_3_13`.

Expected dataset location:

`demos/datasets/reactor_optimization/raw/ex8_3_13`

## Solver

Repository NLP backend.

## Current Backend Status

The current repository solver stack does not yet provide a general nonlinear-programming (NLP) model representation or NLP solver engine.

Therefore, this demo is an explicit backend-boundary demonstration.

It does not:

- silently linearize the nonlinear model,
- convert the problem into a QP,
- convert the problem into an LP,
- or report a solution from a different problem class.

## Workflow

When an NLP backend is available:

1. Load the MINLPLib instance.
2. Construct the nonlinear model.
3. Validate variables, constraints, and nonlinear expressions.
4. Solve using the NLP engine.
5. Report objective value, feasibility, and solver status.

Until then, the demo reports the unsupported NLP backend explicitly.
