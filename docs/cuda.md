# CUDA Backend

Optional GPU execution for the two engines whose iterations are dominated by
sparse products and vector updates: **PDLP** (the whole PDHG iteration runs on
the device) and **QP/ADMM** (a *hybrid*: vector work and sparse products on the
device, the KKT factorisation and solves on the CPU).

> **Status: built, tested, sanitised and benchmarked on one NVIDIA
> configuration** (Windows, MSVC, CUDA 13.4, RTX 5050 Laptop GPU; see
> [Verification status](#verification-status)). Other GPUs, operating systems and
> toolkits are expected to work but have not been tested. The `Auto` thresholds
> are calibrated on that one machine and are heuristics, not guarantees.

---

## 1. Scope

| Engine | On the GPU | On the CPU |
| :--- | :--- | :--- |
| PDLP | Trial step (both fused halves), linesearch reductions, commit, iterate averaging, restart-from-average, activity refresh | Ruiz scaling, preconditioner, spectral-norm estimate (once, before iterating); linesearch/restart/step-size decisions (scalars); termination, infeasibility certificates and polishing *scoring* at check boundaries |
| QP (ADMM) | `A x`, `Aᵀ v`, `P x` (cuSPARSE), right-hand-side assembly, box projection + dual ascent, residual norms and objective | KKT factorisation and triangular solves, rho adaptation, termination test, certificates, polishing |
| Dual simplex, barrier, branch-and-cut (MILP), MIQP, NLP (elastic SQP) | — | Everything. These engines have no CUDA backend; a CUDA request runs them on the CPU and says so (section 4). |

Only PDLP and ADMM QP have CUDA code. Branch-and-cut solves its node LPs with
the dual simplex. MIQP solves its node QPs with ADMM under default options
(`Auto`, whose QP threshold never selects CUDA), so MILP and MIQP solves run
entirely on the CPU.

Presolve, classification, dispatch and postsolve are untouched. A CUDA backend is
a compute device for one engine call on the already-presolved model; it never
presolves, rescales outside the engine, or changes what the engine returns.

## 2. Prerequisites

* NVIDIA GPU with a driver supporting your toolkit.
* CUDA Toolkit providing `nvcc`, the CUDA runtime and cuSPARSE.
* CMake ≥ 3.20 (the project minimum). `CMAKE_CUDA_ARCHITECTURES=native` as a
  default needs CMake ≥ 3.24; with older CMake, pass an explicit architecture list.
* **Windows:** `nvcc` only works with the MSVC host compiler, so the whole
  project must be configured with the Visual Studio toolchain. MinGW/MSYS2 g++
  builds cannot enable CUDA.

### Tested configuration

The only configuration the CUDA path has been built and run on:

| | |
| :--- | :--- |
| OS | Windows 11 |
| Host compiler | Visual Studio 2022 Build Tools (MSVC), Visual Studio CMake generator |
| CUDA Toolkit | 13.4 (`nvcc` V13.4.92) |
| GPU | NVIDIA GeForce RTX 5050 Laptop GPU, compute capability 12.0, 8 GB, driver 610.74 |
| Architecture | `-DCMAKE_CUDA_ARCHITECTURES=120` (`sm_120`) |
| Precision | `double` throughout |

### Minimum toolkit (inferred, not tested)

Derived from the APIs the code uses, not from a tested build:

* C++17 device code (`CUDA_STANDARD 17`): CUDA 11.0 or newer.
* Generic cuSPARSE SpMV with `CUSPARSE_SPMV_CSR_ALG2` (deterministic CSR
  product): CUDA 11.2 or newer.
* A GPU's architecture must be supported by the toolkit; e.g. Blackwell
  (compute capability 12.x, such as the RTX 50 series) requires CUDA 12.8+.

Libraries linked: `CUDA::cudart_static` (both engines) and `CUDA::cusparse`
(QP only). cuBLAS and cuSOLVER are not used.

## 3. Building

CUDA is **off by default**. A default configure never looks for a CUDA compiler.

```bash
# Whole project, both engines
cmake -S . -B build-cuda -DCMAKE_BUILD_TYPE=Release -DOPTIMSOLVER_ENABLE_CUDA=ON
cmake --build build-cuda -j
ctest --test-dir build-cuda --output-on-failure

# A machine without a GPU (CI, packaging) must name the architectures:
cmake -S . -B build-cuda -DOPTIMSOLVER_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="80;86;89;90"

# Engines standalone (each builds its own tests, including the CUDA ones)
cmake -S pdlp_engine -B build-pdlp-cuda -DPDLP_ENABLE_CUDA=ON
cmake -S qp_engine   -B build-qp-cuda   -DQP_ENABLE_CUDA=ON
```

| Option | Default | Effect |
| :--- | :--- | :--- |
| `OPTIMSOLVER_ENABLE_CUDA` | `OFF` | Enables CUDA at the top level and defaults both engine options below to `ON`. |
| `PDLP_ENABLE_CUDA` | `OPTIMSOLVER_ENABLE_CUDA`, else `OFF` | Builds `pdlp_cuda` and `pdlp_cuda_tests`. |
| `QP_ENABLE_CUDA` | `OPTIMSOLVER_ENABLE_CUDA`, else `OFF` | Builds `qp_cuda` and `qp_cuda_tests`. |
| `OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE` | `OFF` | CUDA tests fail instead of skipping if no usable GPU exists; requires a CUDA backend. |
| `CMAKE_CUDA_ARCHITECTURES` | `native` (CMake ≥ 3.24) | Never hard-coded; any explicit value wins. |

Build-system guarantees:

* With CUDA off, no `.cu` file is compiled and a stub (`src/cuda/cuda_unavailable.cpp`)
  provides the same API, reporting that CUDA was not compiled in.
* Device code lives in separate static targets (`pdlp_cuda`, `qp_cuda`), so the
  engines' host flags (`-march=native`, `-Wpedantic`, LTO) never reach `nvcc`,
  and no CUDA flag reaches ordinary C++ targets.
* IEEE semantics are enforced: `--prec-div=true --prec-sqrt=true --ftz=false`,
  and configuring with `--use_fast_math` (or `ftz=true`, `prec-div=false`,
  `prec-sqrt=false`) in `CMAKE_CUDA_FLAGS` is a configure error.
* FMA contraction follows the host compiler of each engine:
  * **PDLP under MSVC:** `pdlp_cuda` is compiled with `-fmad=false`, because
    MSVC's `/fp:precise` never contracts `a*b+c`. With contraction on, the
    device rounds once where the CPU rounds twice, and on diverging
    (infeasible) PDHG paths that difference grew until the backends took
    different restart and certificate decisions.
  * **PDLP with GCC/Clang:** nvcc's default `-fmad=true` stays on, matching the
    host build's `-ffp-contract=fast`.
  * **QP:** `qp_cuda` keeps nvcc's default `-fmad=true` on every platform. Its
    CPU/CUDA comparisons are tolerance-based (section 8), not bitwise.
* All solver arithmetic stays in `double`.
* Enabling an engine CUDA backend registers its CUDA tests even when
  `PDLP_BUILD_TESTS` or `QP_BUILD_TESTS` is `OFF`. Ordinary engine unit tests
  still follow those switches.

## 4. Selecting a backend

Every layer exposes the same three-valued choice, defaulting to `Auto`:

| Layer | Field / flag |
| :--- | :--- |
| PDLP | `pdlp::PdlpOptions::backend`, `cudaDevice`, `cudaNonzeroThreshold` |
| QP | `qp::AdmmOptions::backend`, `cudaDevice`, `cudaNonzeroThreshold` |
| Pipeline | `solver::SolverOptions::backend`, `cudaDevice` |
| CLI | `--backend auto\|cpu\|cuda`, `--cuda-device <index>` |

| Value | Behaviour |
| :--- | :--- |
| `cpu` | Always the CPU implementation. |
| `cuda` | The CUDA backend or an error. If the build lacks CUDA, no device is usable, the device index is out of range, or setup fails (e.g. insufficient device memory), the solve **fails with the reason** -- `InvalidProblem` from the engine, `Unsupported` from the pipeline, a non-zero exit from the CLI. It is never silently run on the CPU. |
| `auto` | CUDA only when all hold: the build has CUDA, the problem is at least `cudaNonzeroThreshold` nonzeros, the device is usable, and setup succeeds. Otherwise the CPU, with the reason recorded. |

Selection is deterministic: it depends only on the options, the build, the
device, and the size of the (scaled) problem.

What actually ran is always recorded, never inferred from the request:

* `PdlpResult::executedBackend` / `backendMessage`, `AdmmResult::executedBackend` / `backendMessage`
* `SolveResult::executedBackend` / `backendReason`
* The CLI prints `Compute backend: …` when `--backend` was given or a GPU ran
  (the default output is unchanged).

Engines without a CUDA backend (dual simplex, barrier, branch-and-cut, MIQP,
NLP) run on the CPU whatever was requested; with `--backend cuda` the result
says so (`"dual_simplex has no CUDA backend; ran on the CPU"`).

### CPU fallback, case by case

| Situation | `cpu` | `auto` | `cuda` |
| :--- | :--- | :--- | :--- |
| Build without CUDA (`OPTIMSOLVER_ENABLE_CUDA=OFF`) | CPU | CPU; reason names the missing CUDA support | Refused: `Unsupported` (pipeline), `InvalidProblem` (engine), non-zero CLI exit |
| CUDA build, no usable device / bad device index | CPU | CPU; reason is the device probe's message | Refused, as above |
| Problem below `cudaNonzeroThreshold` | CPU | CPU; reason quotes the nonzero count and threshold | CUDA (the threshold applies only to `auto`) |
| Device setup fails (for example, not enough device memory) | CPU | CPU; reason quotes the setup error | Refused, as above |
| Engine has no CUDA backend | CPU | CPU | CPU; reason says the engine has no CUDA backend |
| CUDA error during the solve | — | `NumericalFailure` with the CUDA error; no CPU re-run | same |

`auto` tests the size threshold before it probes the device, so a small problem
reports the threshold as its reason even in a build or machine without CUDA.
Nothing ever falls back silently: the executed backend and reason are part of
every result.

JSON reports (including NLP reports) contain a top-level `compute_backend`:

```json
"compute_backend": {
  "requested": "cuda",
  "executed": "cuda",
  "requested_device": 0,
  "executed_device": 0,
  "reason": "explicit CUDA request"
}
```

`executed` and `executed_device` are `null` when no backend ran (for example,
setup refusal). CPU execution has a null `executed_device`, even when CUDA was
requested. The reason records selection, fallback, or refusal; it is not an
inference from the requested option.

### Auto thresholds

| Engine | Default `cudaNonzeroThreshold` | Basis |
| :--- | :--- | :--- |
| PDLP | 1,000,000 | Measured (section 9). Against a well-configured CPU run (12 threads), CUDA is slower up to 500k nonzeros, first faster at 750k (1.17×), and 1M is the first size with a clear margin (1.49×), widening to 2.9× at 5M. Kept unchanged. |
| QP | `INT64_MAX` (Auto never picks CUDA) | Measured (section 9): the hybrid backend ranged from 0.45× to 1.02× of the CPU from 19k to 283k nonzeros, because the KKT factorisation and solves stay on the CPU and dominate. `backend = cuda` always runs it. |

The thresholds are heuristics calibrated on one machine (section 2), not
performance guarantees. A machine with a stronger CPU, a weaker GPU or a
different cache size will have a different crossover. Re-run the benchmarks in
section 9 before changing them.

## 5. PDLP on the device

**Residency.** The scaled problem -- CSR *and* CSC of `A` (64-bit offsets,
32-bit indices, never truncated), bounds, objective, both diagonal
preconditioners -- is uploaded once when the backend is created. The iterate
`x, y`, the carried activity `A x`, the three trial buffers and the running
average stay on the device for the whole solve.

**The step.** The CPU kernel's fused structure is kept:

* primal half over CSC columns: `(Aᵀy)_j → gradient → prox/clip → x'_j`, plus
  `Σ dx²/T`;
* dual half over CSR rows: `(A x')_i → extrapolation → Moreau prox → y'_i, (A x')_i`,
  plus `Σ dy²/Σ` and `Σ dy·(A dx)`.

The kernel boundary between the halves is the one global dependency PDHG has.
The per-coordinate arithmetic is **the same source** on both backends
(`pdlp/pdhg_math.h`, `__host__ __device__`), including `std::min/std::max`
comparison order: it differs from `fmin/fmax` on NaN, and NaN propagation is
how numerical failure is detected.

**Work mapping.** A line (row or column) is reduced by a group of G lanes,
G ∈ {1,2,4,8,16,32} chosen from the matrix's mean line length (largest power of
two not above it). Lines with more than `heavyLineNonzeros` (default 4096)
nonzeros are excluded from the group kernel and reduced by one thread block
each, so a few dense rows or columns cannot stall a warp. Both are heuristics
awaiting profiling (`CudaBackendConfig` overrides them).

**Reductions.** No floating-point atomics. Each block writes one partial per
sum; a single-block pass adds them in index order. Results are therefore
reproducible run-to-run on a given device and launch configuration (unlike the
multi-threaded CPU path, whose dynamic chunk claiming makes its reductions -- and
occasionally its iteration count -- vary between runs).

**Per-iteration traffic.** One host synchronisation per linesearch trial, to
read three doubles (`Σdx²/T`, `Σdy²/Σ`, the interaction). The host needs them to
accept or reject the step. Commit is a buffer swap; averaging is a device
kernel with the blend fraction computed on the host in the same order as the
CPU, so the fractions are bitwise identical.

**Check boundaries (every `terminationCheckFrequency` iterations, default 100).**
The host downloads `x`, `y`, and the average (at most `2(n+m)` doubles each) and
runs the *existing* CPU termination checker, infeasibility detector and restart
policy on them, multi-threaded. This is the correctness-first stage: the checks
are the identical code the CPU path runs, not a port. Its cost is reported as
`PdlpResult::hostCheckSeconds` (of which `backendProfile.snapshotSeconds` is the
download) so a benchmark can decide whether GPU-side termination is worth
building. Restart-from-average is a device-to-device copy plus an activity
refresh; the host-side restart anchor is taken from the snapshot already on the
host.

**Polishing** reuses the same device backend (the problem is already resident).

**Failure handling.** A CUDA error mid-solve returns `NumericalFailure` with
`statusMessage = "CUDA backend error: <API, file:line, status>"`; device memory
is released by RAII. The up-front memory check refuses a problem that would not
fit (`Auto` then uses the CPU; `Cuda` reports it).

## 6. QP on the device (hybrid)

| | Where | Per iteration |
| :--- | :--- | :--- |
| `A`, `Aᵀ` (as its own CSR), `P` (full symmetric CSR, as on the host) | device, uploaded once | — |
| `z, y, A x, z_old, y_old, best y` | device | — |
| `x` | host (KKT output) and device copy | upload `n` doubles |
| x-update right-hand side | built on device | download `n` doubles |
| KKT factorisation + triangular solves | **CPU** (`KktSolver`, unchanged) | — |
| `‖Ax−z‖`, `‖ρAᵀΔz‖` (or `‖Px+q‖`), objective | device reduction | download 3 doubles |
| `y, y_old, z, A x` for the termination/certificate checks | device | download at check boundaries only |

Every sparse product is a non-transposed, deterministic CSR SpMV
(`CUSPARSE_SPMV_CSR_ALG2`); `Aᵀ` is stored explicitly because a transposed CSR
product in cuSPARSE accumulates with atomics. Descriptors and SpMV workspaces
are created once. Offsets and indices are stored 32-bit when `nnz < 2³¹`
(provably safe), otherwise 64-bit; offsets are never truncated.

Temporary CSR index conversion buffers use `uploadNewAndWait`: their stream
finishes the upload before the host vector is destroyed. This adds setup-only
synchronisations, counted in the backend profile, and does not depend on
pageable-memory staging behaviour of `cudaMemcpyAsync`.

Rho adaptation, refactorisation, the termination test, both certificates and
polishing run on the host exactly as before. `AdmmResult::kktSolveSeconds` /
`kktFactorSeconds` (measured on every backend) and `backendProfile` quantify the
split. **Expect the KKT solve to bound any speed-up** (Amdahl): the GPU can only
remove the sparse-product and vector share of each iteration. Moving the KKT
solve to the GPU (cuSOLVER dense Cholesky for small dense KKT systems, or a GPU
sparse factorisation) is deliberately not attempted until profiling justifies it.

## 7. Known limitations

* Tested on a single configuration (section 2). Other GPUs, Linux and other
  toolkits are untested.
* PDLP termination, certificates and restart scoring run on the host at check
  boundaries (correctness-first stage; see section 5).
* Ruiz scaling, the diagonal preconditioner and (without preconditioning) the
  power iteration run on the CPU before upload.
* Group size, heavy-line threshold and grid sizing are unprofiled heuristics.
  The `Auto` thresholds come from one machine's benchmarks (section 9).
* One GPU per solve (`cudaDevice`); no multi-GPU partitioning.
* QP is hybrid; its KKT solve is on the CPU and `Auto` never selects it. In
  the measured range it was never meaningfully faster than the CPU.
* PDLP's CUDA path still runs termination checks on the host with
  `threadCount` threads, so it inherits the CPU threading behaviour described
  in section 9 at check boundaries.
* On Windows, CUDA builds require MSVC for the whole project; the MinGW build
  this repository is usually developed with cannot enable CUDA.

## 8. Tests

| Test | Builds | What it checks |
| :--- | :--- | :--- |
| `pdlp_tests` (new cases) | `PDLP_BUILD_TESTS=ON` | shared math NaN/∞ semantics; backend contract CPU-vs-CPU on 10 matrix shapes; explicit CUDA never falls back; Auto resolution; forced CPU is bitwise the default; invalid device rejected |
| `qp_tests` (new cases) | `QP_BUILD_TESTS=ON` | ADMM backend checks; independent original-model KKT checks across fixture/option variants (solve tolerance 1e-8, componentwise check 1e-5); deliberate corrupted-result rejection |
| `test_compute_backend` | always | pipeline plumbing, `Unsupported` for an unavailable CUDA request, reasons recorded, CUDA requests on dual simplex and barrier acknowledged as CPU runs |
| `test_cli` (new cases) | always | `--backend` / `--cuda-device` parsing, refusal, reporting |
| `backend_json_provenance` | Python available | parsed CLI JSON for CPU/Auto selection, explicit CUDA refusal, and CPU-only LP/NLP engines |
| `pdlp_cuda_tests` | CUDA only | CPU-vs-CUDA kernel contract for every group size, with and without forced heavy lines, on 10 shapes (empty, one-nonzero, `m = 0`, tall, wide, skewed, 1e±9 scaled, all bound kinds); rejected-trial correctness; averaging; restart; NaN detection; full solves over 8 option variants × feasible/infeasible/unbounded; bitwise run-to-run reproducibility of CUDA solves (status, iteration and trial counts, objective, solution, duals and rays); Auto selection; bad device index |
| `qp_cuda_tests` | CUDA only | CPU-vs-hybrid backend contract (lock-stepped on the same KKT solutions, including a rho change); full solves over 5 variants × 7 models; independent original-model primal feasibility, stationarity, multiplier sign, complementarity and objective checks for optimal results; temporary 32/64-bit upload lifetime; bitwise run-to-run reproducibility of hybrid solves on all 7 models |

The CUDA tests exit with code 77, which CTest reports as **skipped**, when no
usable device exists, so a CUDA build on a GPU-less machine is not a failure.
For a GPU validation job, configure `OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE=ON`;
missing hardware then fails the tests. Run `ctest -L cuda` to select both GPU
suites. Direct executable invocation can use the environment variable
`OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE=1` for the same strict behaviour.

## 9. Benchmarks

### Harnesses

* `pdlp_bench [rows] [cols] [nnzPerRow] [iterationLimit] [threads] [cpu|cuda|auto] [repeats]`
  generates a feasible LP from a fixed seed (12345), performs one untimed
  warm-up solve and `repeats` timed solves, and reports the median: end-to-end
  time (everything inside `PdlpSolver::solve`), device setup/upload, host
  check and snapshot time, bytes each way, host synchronisations and
  µs/iteration. The last output line is `json {...}` with the same figures and
  every timed run.
* `qp_bench sweep <cpu|cuda|auto> <n> [nnzPerRow] [repeats] [threads]` does the
  same for a generated banded convex QP (seed 2024): `n` variables, `n`
  constraint rows with `nnzPerRow` nonzeros in a band around the diagonal, a
  diagonally dominant tridiagonal `P`, and rows that cycle through equalities,
  one-sided and ranged constraints. It times `QpSolver::solve` end to end and
  adds the KKT solve/factor split. `qp_bench [cpu|cuda|auto]` without `sweep`
  runs the original small fixed cases.
* `benchmarks/cuda/run_cuda_benchmarks.ps1 -BuildDir <cuda build>` runs both
  sweeps, CPU and CUDA back to back per instance, and writes `pdlp.json` /
  `qp.json` (every raw record), `pdlp.csv` / `qp.csv` (CPU next to CUDA per
  instance) and `summary.md`. `-Threads N` fixes the solver thread count;
  `-Resume` reuses finished instances after an interruption.

Speedup = CPU median / CUDA median. Objective error = |CPU objective − CUDA
objective|. CUDA times include device setup, host-to-device upload, every
kernel and reduction, host-side checks and the final download; nothing is
excluded.

### Results (2026-10-04)

Machine as in section 2, plus an Intel Core i7-13645HX (6 performance + 8
efficiency cores, 20 hardware threads). Release build of the root project
with MSVC (`/O2 /fp:precise`), double precision, 1 untimed warm-up and 5
timed end-to-end solves per backend and instance, median reported. The
laptop was in ordinary desktop use (browser and editor open), not an isolated benchmark
host. Raw data: [`benchmarks/results/cuda/`](../benchmarks/results/cuda/).

**PDLP, CPU fixed at 12 threads.** Generated LPs with `rows × 2·rows` columns
and 10 nonzeros per row, solved to the default tolerances (1e-6) with
polishing off. Every solve on both backends ended `optimal`.

| Instance | NNZ | CPU (s) | CUDA (s) | Speedup | Objective error | Iterations CPU / CUDA | CPU min–max (s) | CUDA min–max (s) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 × 2,000 | 9,977 | 0.074 | 0.359 | 0.21× | 1.6e-05 | 6000 / 6000 | 0.071–0.077 | 0.351–0.368 |
| 10,000 × 20,000 | 99,976 | 0.464 | 1.262 | 0.37× | 2.2e-05 | 11200 / 11200 | 0.461–0.471 | 1.254–1.294 |
| 25,000 × 50,000 | 249,972 | 0.975 | 1.493 | 0.65× | 6.9e-05 | 7800 / 7800 | 0.830–1.065 | 1.464–1.530 |
| 50,000 × 100,000 | 499,974 | 2.313 | 2.529 | 0.91× | 5.5e-05 | 8000 / 8000 | 2.288–2.331 | 2.499–2.539 |
| 75,000 × 150,000 | 749,971 | 4.764 | 4.073 | 1.17× | 2.6e-05 | 9000 / 9000 | 4.711–4.838 | 4.021–4.183 |
| 100,000 × 200,000 | 999,977 | 8.178 | 5.480 | 1.49× | 4.3e-06 | 9000 / 9000 | 8.141–8.497 | 5.389–5.485 |
| 150,000 × 300,000 | 1,499,978 | 16.746 | 7.667 | 2.18× | 1.5e-04 | 8600 / 8800 | 16.058–18.544 | 7.611–7.767 |
| 200,000 × 400,000 | 1,999,974 | 24.762 | 9.845 | 2.52× | 2.6e-05 | 8600 / 8600 | 24.629–24.906 | 9.731–10.005 |
| 500,000 × 1,000,000 | 4,999,978 | 75.848 | 26.256 | 2.89× | 1.7e-04 | 8800 / 8800 | 75.078–79.043 | 25.946–32.023 |

Objective errors are at most 6e-9 relative to the objective (|objective| is
about 2.8e3 to 1.4e6 on these instances), well within the 1e-6 solve tolerance.

**PDLP, CPU with all 20 hardware threads (the default `threadCount = 0`).**
Same instances and method, run first. The CPU timings were bimodal: the same
solve, with the same iteration count, took anywhere from 1.04 s to 8.34 s at
250k nonzeros and from 85 s to 353 s at 5M. A separate scan at 500k nonzeros
(3 runs per setting) measured 11.96–12.42 s with 1 thread, 2.53–3.13 s with
6, 2.75–2.91 s with 8, 2.42–2.48 s with 12, and 2.74 s to 817 s with 14. The
CPU kernel's synchronisation does not tolerate threads being descheduled on
this hybrid-core CPU; that is a property of the CPU engine on this machine,
not of the CUDA path, and it is why the 12-thread series above is the basis
for the threshold. The CUDA medians here are also slower than in the
12-thread series, because the host-side checks use the same thread count.

| Instance | NNZ | CPU (s) | CUDA (s) | Speedup | Objective error | CPU min–max (s) | CUDA min–max (s) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 × 2,000 | 9,977 | 0.081 | 0.454 | 0.18× | 1.6e-05 | 0.079–0.091 | 0.451–0.484 |
| 10,000 × 20,000 | 99,976 | 0.793 | 2.251 | 0.35× | 1.6e-05 | 0.645–7.323 | 1.579–2.452 |
| 25,000 × 50,000 | 249,972 | 8.080 | 2.390 | 3.38× | 8.5e-06 | 1.043–8.343 | 1.684–2.479 |
| 50,000 × 100,000 | 499,974 | 11.274 | 2.781 | 4.05× | 3.6e-05 | 2.224–16.370 | 2.677–3.082 |
| 75,000 × 150,000 | 749,971 | 5.813 | 5.929 | 0.98× | 2.2e-04 | 5.148–29.732 | 5.908–5.951 |
| 100,000 × 200,000 | 999,977 | 40.438 | 7.657 | 5.28× | 1.3e-05 | 39.009–40.956 | 7.133–7.702 |
| 150,000 × 300,000 | 1,499,978 | 59.545 | 10.667 | 5.58× | 2.7e-04 | 51.920–61.052 | 10.498–10.695 |
| 200,000 × 400,000 | 1,999,974 | 80.831 | 17.755 | 4.55× | 5.6e-05 | 63.262–82.612 | 10.657–29.563 |
| 500,000 × 1,000,000 | 4,999,978 | 149.773 | 37.739 | 3.97× | 7.6e-05 | 84.948–352.515 | 37.596–37.927 |

**Threshold decision.** With the stable CPU baseline the crossover lies
between 500k and 750k nonzeros, and CUDA's margin grows steadily above it.
The default `cudaNonzeroThreshold = 1,000,000` is the first measured size where
CUDA wins with a clear margin against both CPU configurations, so it is kept.
Lowering it to about 750k would gain at most 17% at those sizes on this
machine, and a stronger CPU would move the crossover up, not down.

**QP, CPU vs hybrid CUDA, all threads.** Generated banded QPs, `n` variables
and `n` rows with 8 nonzeros each, default ADMM options. Every solve ended
`Optimal`, with the same iteration count on both backends.

| Instance (n) | NNZ (P+A) | CPU (s) | Hybrid CUDA (s) | Speedup | CUDA KKT solve (s) | CUDA setup (s) | H2D / D2H (MB) | Objective error |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2,000 | 18,869 | 0.023 | 0.051 | 0.45× | 0.007 | 0.001 | 2.7 / 2.5 | 0 |
| 5,000 | 47,209 | 0.125 | 0.152 | 0.82× | 0.020 | 0.001 | 6.8 / 6.2 | 4.0e-14 |
| 10,000 | 94,375 | 0.448 | 0.500 | 0.90× | 0.042 | 0.002 | 13.6 / 12.4 | 5.7e-13 |
| 20,000 | 189,035 | 5.514 | 5.771 | 0.96× | 0.749 | 0.008 | 42.5 / 41.2 | 8.0e-13 |
| 30,000 | 283,425 | 11.816 | 11.579 | 1.02× | 1.165 | 0.010 | 63.7 / 61.8 | 1.0e-12 |

The hybrid backend is not faster in this range. Most of the time is the CPU
KKT factorisation, the same on both backends and growing faster than linearly
with `n` for this sparse Cholesky: 0.08 s, 0.36 s, 4.1 s and 9.6 s of the
CPU solve at n = 5,000, 10,000, 20,000 and 30,000 (`kkt_factor_seconds` in
`qp.json`). Moving the remaining work to the GPU cannot reduce that, so `Auto`
keeps QP on the CPU.

### CPU baseline on the development machine (earlier, MinGW)

Intel Core i7-13645HX (14 cores / 20 threads), MSYS2 UCRT64 g++ 15.2, standalone
Release build (`-O3 -march=native`, LTO), `iterationLimit = 1000`, polishing off,
3 timed runs after 1 warm-up. The laptop throttles under sustained load:
warm-up and median can differ by 2–3× on the largest size, so both are shown.

| rows × cols | nnz | threads | warm-up (s) | median e2e (s) | µs / iteration (median) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| 2,000 × 4,000 | 15,983 | 20 | 0.043 | 0.035 | 34.7 |
| 2,000 × 4,000 | 15,983 | 1 | 0.036 | 0.034 | 33.8 |
| 20,000 × 40,000 | 199,978 | 20 | 0.141 | 0.132 | 131.6 |
| 20,000 × 40,000 | 199,978 | 1 | 0.749 | 0.765 | 765.3 |
| 100,000 × 200,000 | 999,977 | 20 | 1.140 | 3.263 | 3,263.3 |
| 100,000 × 200,000 | 999,977 | 1 | 11.507 | 20.367 | 20,366.8 |

The smallest case is below `parallelNonzeroThreshold`, so it runs serially
either way. These figures come from a different compiler and build (MinGW g++,
`-march=native`, LTO, a fixed 1,000 iterations) than the MSVC results above,
so the two sets are not directly comparable.

## 10. Validating on a CUDA machine

The scripts in `tools/cuda/` run every check below and write logs plus a
summary. They are PowerShell (Windows PowerShell 5.1 or PowerShell 7). Each
child process is started directly, so a test that fails is reported
separately from an executable the operating system refused to start (for
example, Windows Application Control blocking a freshly built binary).

```powershell
# Everything: fresh CPU-only and CUDA Release builds, full CTest on both,
# ctest -L cuda, direct runs of both CUDA test binaries, Compute Sanitizer
# (memcheck, initcheck; -Racecheck adds racecheck) and 10 repeated runs.
tools/cuda/validate_cuda.ps1 -CudaArchitectures 120
#   -> build-validate-logs/validation-summary.md

# Individual steps against an existing CUDA build
tools/cuda/run_cuda_sanitizers.ps1     -BuildDir build-validate-cuda
tools/cuda/check_cuda_determinism.ps1  -BuildDir build-validate-cuda -Runs 10

# Benchmarks (needs -DPDLP_BUILD_TOOLS=ON -DQP_BUILD_BENCH=ON; validate_cuda.ps1 sets both)
benchmarks/cuda/run_cuda_benchmarks.ps1 -BuildDir build-validate-cuda -Threads 12
```

`validate_cuda.ps1` sets `PYTHONUTF8=1` for the test runs, so the Python-driven
CLI tests decode output as UTF-8 rather than the Windows ANSI code page. Tests
listed in `-KnownFailures` (by default `nlp_elastic_kkt_tests`, a pre-existing
failure unrelated to CUDA) are reported separately, not hidden.

**Sanitizers.** Every tool must report `ERROR SUMMARY: 0 errors` and the test
binary itself must pass. memcheck runs with `--leak-check full`. racecheck is
opt-in (`-Racecheck`) because it is slow: over the full `pdlp_cuda_tests`
suite it ran for more than 50 minutes on the tested GPU without finishing.
Run it when kernel or synchronisation code changes.

**Determinism.** Both CUDA test binaries assert bitwise run-to-run
reproducibility of CUDA solves inside one process. The determinism script
also runs each binary repeatedly and compares, across processes, every
per-case pass/fail line, the CPU and CUDA iteration counts and the transfer
byte counts. Wall-clock figures are removed before comparing; they are never
compared.

The equivalent manual commands, on Linux or from a Visual Studio developer
prompt:

```bash
cmake -S . -B build-cuda -DCMAKE_BUILD_TYPE=Release -DOPTIMSOLVER_ENABLE_CUDA=ON \
      -DOPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE=ON -DCMAKE_CUDA_ARCHITECTURES=120 \
      -DPDLP_BUILD_TESTS=ON -DQP_BUILD_TESTS=ON -DPDLP_BUILD_TOOLS=ON -DQP_BUILD_BENCH=ON
cmake --build build-cuda --config Release -j
ctest --test-dir build-cuda -C Release --output-on-failure
ctest --test-dir build-cuda -C Release -L cuda

compute-sanitizer --tool memcheck --leak-check full <path>/pdlp_cuda_tests
compute-sanitizer --tool initcheck  <path>/pdlp_cuda_tests
compute-sanitizer --tool racecheck  <path>/pdlp_cuda_tests
compute-sanitizer --tool memcheck --leak-check full <path>/qp_cuda_tests
compute-sanitizer --tool initcheck  <path>/qp_cuda_tests

<path>/pdlp_bench 100000 200000 10 200000 12 cpu 5
<path>/pdlp_bench 100000 200000 10 200000 12 cuda 5
<path>/qp_bench sweep cuda 20000 8 5

# Profiles
nsys profile -o pdlp <path>/pdlp_bench 100000 200000 10 2000 0 cuda 1
ncu --set full -k regex:"lineKernel|heavyLineKernel" -c 20 \
    <path>/pdlp_bench 100000 200000 10 200 0 cuda 1
```

What to look at first in a profile: whether `hostCheckSeconds` is a material
share of the PDLP solve (then move termination onto the device); achieved
bandwidth and load balance of `lineKernel` across group sizes; and, for QP,
the KKT factorisation and solve against the total.

## Verification status

### CUDA validation on NVIDIA hardware (2026-10-04)

On the configuration in section 2, with `tools/cuda/validate_cuda.ps1
-CudaArchitectures 120` (fresh build trees):

| Check | Result |
| :--- | :--- |
| CPU-only Release build (`OPTIMSOLVER_ENABLE_CUDA=OFF`) | Built; no CUDA compiler configured. Only pre-existing MSVC warnings (C4324 alignment padding in PDLP headers, D9025 `/UNDEBUG` in tests). |
| CUDA Release build (`sm_120`) | Built `pdlp_cuda`, `qp_cuda`, both CUDA test binaries, `optimsolver`, `pdlp_bench`, `qp_bench`. No warnings from CUDA sources. |
| CTest, CPU-only build | 75 of 76 passed; 1 known pre-existing failure (below) |
| CTest, CUDA build | 77 of 78 passed; 1 known pre-existing failure (below) |
| `ctest -L cuda` | 2 of 2 passed, none skipped (`OPTIMSOLVER_REQUIRE_CUDA_TEST_DEVICE=ON`) |
| `pdlp_cuda_tests` / `qp_cuda_tests`, run directly | Passed: 178 and 49 cases |
| Compute Sanitizer memcheck (`--leak-check full`) | PDLP and QP: `ERROR SUMMARY: 0 errors`, 0 bytes leaked |
| Compute Sanitizer initcheck | PDLP and QP: `ERROR SUMMARY: 0 errors` |
| Repeated runs | Each CUDA test binary 10 times: 10/10 passed, with identical per-case results, iteration counts and transfer sizes in every run; the in-process bitwise reproducibility cases passed |

racecheck was not re-run: no kernel or synchronisation code changed since the
CUDA backend was first validated (the only device-source change is the removal
of an unused host-side helper), and a full-suite racecheck did not finish
within 50 minutes.

Known, unrelated to CUDA:

* `nlp_elastic_kkt_tests` fails: an inner QP stops at `IterationLimit` with
  primal and dual residuals near 1e9. The same failure was previously
  reproduced on the commit immediately before the CUDA backend was merged.
* On the first CTest pass after a fresh build, Windows Application Control
  refused to start a few newly built test executables (CTest `BAD_COMMAND`,
  "Not Run"; a different set in each build tree). On the next pass all of them
  started and passed. `validate_cuda.ps1` reports such tests as launch-blocked
  and retries them directly; they are never counted as passes or as solver
  failures.
* Python-driven CLI tests need `PYTHONUTF8=1` on Windows to decode output
  under the cp1252 code page; `validate_cuda.ps1` sets it.

### Earlier, CPU-only verification

Before CUDA hardware was available, the backend refactor was verified CPU-only:

* **Bitwise-unchanged CPU results:** a fingerprint harness solving 96 serial
  LP/QP configurations (every option toggle, feasible/infeasible/unbounded,
  iteration limits) produced byte-identical output -- statuses, iteration and
  trial counts, objectives, residuals, step sizes, and hashes of every solution
  and ray vector -- before and after the backend refactor.
* The review fixes passed 70/70 CTests on macOS without a CUDA toolkit,
  including the JSON provenance test, and the QP KKT checker rejects
  deliberately corrupted primal, dual, complementarity and objective data.

## Addendum (2026-10-02): CUDA in KAIRO Desktop

*Added for KAIRO v1 (PR #20).*

KAIRO Desktop exposes the same backend selection as the CLI: **Auto / CPU / CUDA** and the CUDA device. It offers these only for engines that have a CUDA backend (automatic dispatch, PDLP and ADMM QP), and passes them unchanged into `SolverOptions`. It makes no backend decision of its own.

What happens follows the record:
- `compute_backend.executed` says what ran.
- `compute_backend.reason` is shown word for word.
- An explicit CUDA request that cannot run is shown as **"Requested backend unavailable"**, not as a solver failure. As section 4 says, such a request is never silently run on the CPU.

The desktop builds produced by the `desktop` GitHub Actions workflow are configured without `OPTIMSOLVER_ENABLE_CUDA`, so they are CPU-only. The desktop's own CUDA path has not been exercised on NVIDIA hardware; the validation above covers the engines, pipeline and CLI underneath it.
