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
- **Pass Rate:** **60 / 60 automated CTest targets pass (100%)**.
- **Assertion Integrity:** Test binaries compile with `-UNDEBUG`, preserving `assert()` statements across both Debug and Release builds.
- **Coverage:** Presolve cascades, dual postsolve reconstruction (29 parameterized test cases), QP matrix conventions, ground-truth MILP enumeration ($2^n$ subsets), CLI options, and process management.

---

## Phase 2 — Interactive Terminal Solver (~60s)

Demonstrate the user-friendly terminal interface for model exploration and interactive inspection.

```bash
./build/optimsolver
```

### Presenter Walkthrough
1. **Startup Banner:** Point out the ASCII mascot banner and supported problem families (`LP · QP · MILP`).
2. **Load Model:** Select Option `[1]` (`Open MPS Model`) and enter:
   ```text
   tests/cli/simple_lp.mps
   ```
3. **Inspect Model Card:** Select Option `[3]` (`Model Information`) to show problem dimensions:
   - 2 variables (both continuous)
   - 1 constraint row
   - Linear objective coefficients
   - Detected problem family: `LP`
4. **Solve Model:** Select Option `[2]` (`Solve Current Model`).
5. **Review Solve Result Card:**
   - Status: `✓ Optimal`
   - Objective value: `8`
   - Iterations: `1`
   - Solved in `< 1 ms`
   - Dual availability: `1 shadow prices, 2 reduced costs`
6. **Exit:** Enter `6` or `q` to return to the shell.

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
```bash
python3 benchmarks/bench.py benchmarks/instances/netlib/afiro.mps \
    --solvers dual_simplex,pdlp,highs --timeout 60 \
    --best-known benchmarks/instances/netlib/best_known.json \
    --out /tmp/afiro_bench.json
```

### Evaluator Evidence
- **Process Isolation:** The benchmark harness (`bench_runner`) executes each solver in its own process group with a POSIX watchdog timer and OS kernel `wait4(..., ru_maxrss)` memory tracking.
- **Parse Cross-Check:** The C++ reader and independent Python reader (`mps_model.py`) verify model equality field-by-field before solving.
- **Independent Checker:** The solution vector is verified by `verify.py` against the **original unscaled model** using frozen tolerances ($10^{-6}$ absolute, $10^{-8}$ relative).
- **Separation of Status and Verdict:**
  - Solver Termination Status: `optimal`
  - Independent Checker Verdict: `optimal_verified` (primal feasible, dual sign valid, KKT stationarity satisfied, duality gap $\le 10^{-6}$).

### Precomputed Portfolio Evidence (Fallback)
If running Python scripts live is constrained by time or environment dependencies:
- Point to committed benchmark results in `benchmarks/results/netlib_smoke.json`.
- Highlight the **8 / 8 Netlib LP smoke suite coverage** achieved by the combined internal portfolio:
  - Dual Simplex verifies 6/8 (stalling on dense cycling models `blend` and `share2b`).
  - PDLP verifies 6/8 (stalling at the tight tolerance boundary on `adlittle` and `sc50b`).
  - Together, the two internal engines complement each other to verify 100% of the suite.

---

## Phase 6 — Closing & Sovereign Value (~30s)

Conclude the live presentation with three concise takeaways:

1. **Sovereign & Modular Architecture:** Clean C++17 design that decouples ingestion, presolve, numerical dispatch, and postsolve. Completely independent of closed-source commercial libraries.
2. **Empirical Verification Rigor:** 60/60 automated CTest targets, 8/8 Netlib smoke suite independently verified, and analytical KKT verification for convex QP.
3. **Engineering Integrity & Roadmap:**
   - Current release provides CPU-based LP, MILP, and QP engines.
   - GPU/CUDA acceleration for linear algebra kernels (matrix-vector multiplication in PDLP and ADMM QP) is designed as a modular plug-in and represents the active engineering roadmap.

---

## Live Demo Safety & Fallback Matrix

| Action / Command | Risk Level | Expected Time | Safety Guidance & Fallbacks |
| :--- | :---: | :---: | :--- |
| `ctest --test-dir build --output-on-failure` | **Zero Risk** | ~2–4 seconds | Safe. Guaranteed 60/60 pass across all targets. |
| `./build/optimsolver` (interactive) | **Zero Risk** | User-paced | Safe. Use `tests/cli/simple_lp.mps`. Avoid typing invalid paths. |
| `./build/optimsolver solve tests/cli/simple_lp.mps` | **Zero Risk** | < 5 ms | Instantaneous batch solve. |
| `./build/optimsolver solve tests/cli/presolve_reduction.mps` | **Zero Risk** | < 5 ms | Demonstrates dimension reduction (3 vars $\to$ 2 vars $\to$ 3 vars). |
| `benchmarks/bench.py` on `afiro.mps` | **Low Risk** | ~2 seconds | Requires Python 3 + SciPy for HiGHS reference. If Python is absent, open `benchmarks/results/netlib_afiro.json`. |
| Full Netlib Smoke Suite (`benchmarks/bench.py`) | **Do Not Run Live** | ~15–30 seconds | Do NOT run all 8 instances live during a 5-minute pitch. Point directly to `benchmarks/results/netlib_smoke.json`. |

### What NOT to Claim During the Demo
- ❌ Do NOT claim PDLP solves all Netlib models (it solves 6/8; Dual Simplex solves the remaining 2).
- ❌ Do NOT claim GPU acceleration is currently executing (it is a roadmap architecture item).
- ❌ Do NOT claim MILP solutions have complete independent dual-bound certificates (primal feasibility and integrality are verified).
- ❌ Do NOT claim our solver universally outperforms commercial tools (focus on sovereignty, modularity, and memory efficiency).
