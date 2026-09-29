# SIH26119 — Live Demonstration

---

### Project & Problem Statement Information

- **Organization:** Mangalore Refinery and Petrochemicals Limited (MRPL)
- **Problem Statement ID:** SIH26119
- **Problem Title:** Indigenous GPU-Accelerated Optimization Solver (Sovereign Alternative to Express / CPLEX)
- **Category:** Software (Smart Automation)
- **Target Presentation Window:** 4–5 Minutes

### Video Recording Links

*Upload the recording to YouTube (Unlisted or Public) or Google Drive, then update the links below:*

- **Primary Demonstration Video:**
  *(To be replaced with final recording link upon upload)*
- **Presentation Slide Deck Backup (Google Drive):**
  https://drive.google.com/file/d/1o24vzDsYfrNa9i11FnNQFsK9G2Me1p-o/view

---

## Demo Objective

This live demonstration presents `optimisation_solver`: an open-source, modular C++17 optimization engine developed for industrial planning, scheduling, and resource allocation. The demonstration walks evaluators through standard MPS model ingestion, structural validation, automatic engine dispatch, numerical solution, invertible presolve/postsolve coordinate and dual reconstruction, process-isolated execution, and independent mathematical verification against unscaled original models.

---

## Prerequisites & Environment Setup

Build and run directly from the local repository build tree:

```bash
# 1. Compile the solver executable and test suite in Release mode
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# 2. Verify binary existence
./build/optimsolver --help
```

- **Executable Location:** `./build/optimsolver`
- **Supported Platforms:** Linux / macOS (C++17 compiler, CMake 3.20+)
- **Verification:** All commands below run out-of-the-box from the repository root.

---

## Phase 1 — Automated Test Suite (~30s)

Demonstrate foundational correctness across all pipeline stages before launching interactive solves.

```bash
ctest --test-dir build --output-on-failure
```

### Expected Evaluator Evidence
- **Pass Rate:** **64 / 64 automated CTest targets pass (100%)**.
- **Assertion Integrity:** Test binaries compile with `-UNDEBUG`, so `assert()` statements stay active in Release builds.
- **Coverage:** Presolve cascades, dual postsolve reconstruction (29 parameterized test cases), QP matrix conventions, ground-truth MILP enumeration ($2^n$ subsets), CLI options, process management, and the NLP engine (4 targets).

---

## Phase 2 — Interactive Terminal Solver (~60s)

Demonstrate the user-friendly terminal interface for model exploration and interactive inspection.

```bash
./build/optimsolver
```

### Presenter Walkthrough
1. **Startup Banner:** Point out the ASCII mascot banner and supported problem families (`LP · QP · MILP`).
2. **Load Model:** Select Option `[1]` (`Open MPS Model or NLP Model`) and enter:
   ```text
   tests/cli/simple_lp.mps
   ```
   Once a model is loaded, the menu changes to: `[1]` Solve Current Model, `[2]` Open Another Model, `[3]` Model Information, `[4]` Solver Settings, `[5]` Help, `[6]` Exit.
3. **Inspect Model Card:** Select Option `[3]` (`Model Information`) to show problem dimensions, then press Enter to return:
   - 2 variables (both continuous)
   - 1 constraint row
   - Objective: Minimize, no quadratic terms
   - Detected problem family: `LP`
4. **Solve Model:** Select Option `[1]` (`Solve Current Model`).
5. **Review Solve Result Card:**
   - Status: `✓ OPTIMAL`
   - Engine: `dual_simplex`
   - Objective value: `8`
   - Primal feasibility `✓`; Dual feasibility `✓ (1 shadow prices, 2 reduced costs)`
   - Iterations: `1`
   - Solved in `< 1 ms`
   - Press Enter at the export prompt to skip saving.
6. **Exit:** Enter `6` (with a model loaded) or `q` to return to the shell.

---

## Phase 3 — Batch CLI & Engine Dispatch (~45s)

Demonstrate non-interactive command-line execution for automated workflows and scripting.

### 1. Default Automatic Dispatch
```bash
./build/optimsolver solve tests/cli/simple_lp.mps
```
- The classifier automatically recognizes a continuous LP and routes to `dual_simplex`.
- Output displays structured terminal cards showing model statistics and solve summary.

### 2. Explicit Numerical Engine Override
Demonstrate switching from simplex to first-order Primal-Dual Hybrid Gradient (PDHG):
```bash
./build/optimsolver solve tests/cli/simple_lp.mps --solver pdlp
```
- Demonstrates engine decoupling: the orchestrator routes the identical model to the `pdlp_engine`.
- Both engines reach the identical optimum (`Objective 8`).

### 3. Structured JSON Export for Automated Pipelines
```bash
./build/optimsolver solve tests/cli/simple_lp.mps --json /tmp/solve_record.json
cat /tmp/solve_record.json
```
- Emits machine-readable solve metadata: solver status, wall-clock time, iteration counts, objective values, primal solution vector, constraint duals, and reduced costs.
- Enforces strict null handling: uncomputed fields are `null`, never fabricated zeros.

---

## Phase 4 — Invertible Presolve & Dual Postsolve (~45s)

Demonstrate the **Single-Presolve Invariant**: presolve reduces the problem once, logs transformation metadata (`PresolveResult`), the solver operates in reduced space, and postsolve restores original coordinates and shadow prices.

```bash
./build/optimsolver solve tests/cli/presolve_reduction.mps --output /tmp/solution.txt
```

### Evaluator Evidence
1. **Dimension Reduction Displayed on Terminal:**
   ```text
   Model      3 variables · 1 constraint
   Reduced    2 variables · 1 constraint
   Engine     dual_simplex
   ```
   - Variable `X3` is fixed ($lb = ub = 5$). Presolve logs `X3 = 5` into the audit log and eliminates it.
   - The numerical engine solves a reduced 2-variable problem.
2. **Original Coordinate & Dual Restoration:**
   Inspect the exported solution:
   ```bash
   cat /tmp/solution.txt
   ```
   ```text
   # Solution for PRESOLVE_RED
   # Status: optimal
   # Objective: 28
   X1 4
   X2 0
   X3 5
   # Dual C1 2
   # Reduced cost X1 0
   # Reduced cost X2 1
   # Reduced cost X3 2
   ```
   - All 3 original variables are restored (`X1=4, X2=0, X3=5`).
   - The objective is re-evaluated in original coordinates: $3(4) + 2(0) + 2(5) + 6 = 28$.
   - Constraint duals and reduced costs are mapped back to original constraints and variables.

---

## Phase 5 — Independent Benchmark Verification (~60s)

Demonstrate that the solver is verified against an **independent external checker** rather than trusting its own self-reported numbers.

### Live Command: Benchmark Single Netlib Model
`afiro.mps` is the only Netlib instance committed to the repository, so this command needs no download. `bench.py` expects a Release build in `build-release/`:
```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release && cmake --build build-release -j
python3 benchmarks/bench.py benchmarks/instances/netlib/afiro.mps \
    --solvers dual_simplex,pdlp,highs --timeout 30 --threads 1 \
    --best-known benchmarks/instances/netlib/best_known.json \
    --out /tmp/afiro_bench.json
```
- **HiGHS reference:** The committed results use the **native `highs` executable**. If `highs` is not on `PATH`, `bench.py` falls back to SciPy and labels the row `highs-scipy`. Do not quote that row's memory figures, because they include the Python interpreter.

### Evaluator Evidence
- **Process Isolation:** The benchmark harness (`bench_runner`) executes each solver as a separate process tree with a watchdog timer (SIGTERM, then SIGKILL) and records peak memory from the kernel via `wait4(..., ru_maxrss)`.
- **Parse Cross-Check:** The C++ reader and independent Python reader (`mps_model.py`) verify model equality field-by-field before solving.
- **Independent Checker:** The solution vector is verified by `verify.py` against the **original unscaled model** using frozen tolerances ($10^{-6}$ absolute, $10^{-8}$ relative).
- **Separation of Status and Verdict:**
  - Solver Termination Status: `optimal`
  - Independent Checker Verdict: `optimal_verified` (primal feasible, dual sign valid, KKT stationarity satisfied, duality gap $\le 10^{-6}$).

### Precomputed Portfolio Evidence (Fallback)
If running Python scripts live is constrained by time or environment dependencies:
- Point to the committed results: `benchmarks/results/NETLIB_RESULTS.md` (readable summary), `benchmarks/results/netlib_lp.json` (full records) and `benchmarks/results/netlib_lp_benchmark.csv`.
- Frame the evidence as an **8-instance Netlib LP smoke suite** (30-second timeout, one thread), not the full Netlib collection:
  - **Dual Simplex** verifies 6/8. On `blend` and `share2b` it stops at a limit (`limit_reached`), returning a feasible point well short of the published optimum.
  - **PDLP** verifies 6/8. On `adlittle` and `sc50b` PDLP reports `optimal` and its point is feasible and matches the published objective. The independent checker could not close an optimality proof from its duals, so those two are recorded as `feasible`.
  - **Together**, the two internal engines independently verify all 8 instances. 4 instances are verified by both, and each covers the other's 2 unverified instances. This is the union of separate runs.
  - **The native HiGHS reference** verifies 8/8 and is faster on average (0.009 s vs 0.090 s PDLP and 0.209 s Dual Simplex mean wall time).

---

## Phase 6 — Closing & Sovereign Value (~30s)

Conclude the live presentation with three concise takeaways:

1. **Sovereign & Modular Architecture:** Clean C++17 design that decouples ingestion, presolve, numerical dispatch, and postsolve. Completely independent of closed-source commercial libraries.
2. **Empirical Verification Rigor:** 64/64 automated CTest targets. On the 8-instance Netlib LP smoke suite, every instance is independently verified optimal by at least one internal engine. The two convex QP fixtures pass independent KKT verification.
3. **Engineering Integrity & Roadmap:**
   - Current release provides CPU-based LP, MILP, and QP engines.
   - GPU/CUDA acceleration for linear algebra kernels (matrix-vector multiplication in PDLP and ADMM QP) is designed as a modular plug-in and represents the active engineering roadmap.

---

## Live Demo Safety & Fallback Matrix

| Action / Command | Risk Level | Expected Time | Safety Guidance & Fallbacks |
| :--- | :---: | :---: | :--- |
| `ctest --test-dir build --output-on-failure` | **Zero Risk** | ~5–15 seconds (machine-dependent) | Safe. 64/64 targets passed on the current build. Run it once before the demo to confirm on the presentation machine. |
| `./build/optimsolver` (interactive) | **Zero Risk** | User-paced | Safe. Use `tests/cli/simple_lp.mps`. Avoid typing invalid paths. |
| `./build/optimsolver solve tests/cli/simple_lp.mps` | **Zero Risk** | < 5 ms | Instantaneous batch solve. |
| `./build/optimsolver solve tests/cli/presolve_reduction.mps` | **Zero Risk** | < 5 ms | Demonstrates dimension reduction (3 vars $\to$ 2 vars $\to$ 3 vars). |
| `benchmarks/bench.py` on `afiro.mps` | **Low Risk** | ~2 seconds | Requires Python 3, a `build-release/` build, and the native `highs` executable for a comparable HiGHS row (otherwise it falls back to `highs-scipy`). If unavailable, open `benchmarks/results/NETLIB_RESULTS.md`. |
| Full Netlib Smoke Suite (`benchmarks/bench.py`) | **Do Not Run Live** | Needs a download first | 7 of the 8 instances are not committed and must be fetched with `python3 benchmarks/fetch_netlib.py --set smoke`. Point to `benchmarks/results/NETLIB_RESULTS.md` / `netlib_lp.json` instead. |

### What NOT to Claim During the Demo
- ❌ Do NOT claim PDLP is verified on all smoke-suite instances. It is verified on 6/8, and Dual Simplex is verified on the remaining 2.
- ❌ Do NOT claim results on "Netlib" in general. The evidence is an 8-instance smoke suite.
- ❌ Do NOT claim GPU acceleration is currently executing (it is a roadmap architecture item).
- ❌ Do NOT claim MILP solutions have complete independent dual-bound certificates (primal feasibility and integrality are verified).
- ❌ Do NOT claim infeasibility or unboundedness results are independently certified (the checker marks them `not_applicable`).
- ❌ Do NOT claim we are faster than HiGHS. The native HiGHS reference was faster on average on this suite.
- ❌ Do NOT claim universally lower memory. Mean peak RSS was lower on this suite, but Dual Simplex peaked at 9.52 MB on `degen2` against 5.34 MB for HiGHS.
- ❌ Do NOT claim our solver universally outperforms commercial tools (focus on sovereignty, modularity, and independent verification).
