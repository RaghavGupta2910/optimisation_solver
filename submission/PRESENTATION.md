# SIH 2026 Project Presentation Slide Deck

---

### Project & Problem Statement Information

- **Organization:** Mangalore Refinery and Petrochemicals Limited (MRPL)
- **Problem Statement ID:** SIH26119
- **Problem Statement Title:** Indigenous GPU-Accelerated Optimization Solver (Sovereign Alternative to Express / CPLEX)
- **Category:** Software
- **Theme:** Smart Automation
- **Repository:** `https://github.com/RaghavGupta2910/optimisation_solver`

### Slide Deck Cloud Links

- **Slide Deck (Google Drive Backup):**
  https://drive.google.com/file/d/1o24vzDsYfrNa9i11FnNQFsK9G2Me1p-o/view

---

## Slide Structure & Technical Presentation Content

### Slide 1: Title & Team Overview
- **Title:** `optimsolver` — Indigenous Modular Mathematical Optimization Solver
- **SIH Context:** Problem Statement SIH26119, sponsored by Mangalore Refinery and Petrochemicals Limited (MRPL).
- **Team Leadership & Engineering Roles:**
  - **Raghav Gupta** — Team Lead — System Architecture, Solver Orchestration & CLI
  - **Harehar Narayan Seth** — Lead Numerical Solver Development & Optimization Algorithms
  - **Aryan Kumar** — Solver Engine Development & Numerical Methods
  - **Preetish Attray** — MILP / Branch-and-Cut Development
  - **Kabir Pahwa** — Postsolve, Dual Reconstruction & Solution Validation
  - **Sarisha Jhinghan** — Presolve, Model Processing & Testing
- *Visual:* MRPL & SIH partner branding, team grid, and repository badges (C++17, CMake 3.20+, CTest 64/64, MIT).

---

### Slide 2: Problem Context & Sovereign Imperative
- **Industrial Workloads:** Refinery operations depend on continuous mathematical optimization for crude distillation planning, fluid catalytic cracking blending, pipeline scheduling, and resource allocation.
- **Sovereign Alternative Motivation:** Commercial solvers (FICO Xpress, IBM CPLEX) are proprietary, closed-source, and impose steep recurrent licensing fees and foreign technology dependencies.
- **Architectural Goal:** Develop an indigenous, modular C++17 optimization solver founded on transparent mathematical contracts, clean subsystem separation, and reproducible verification.
- *Visual:* Diagram contrasting proprietary monolithic solver stacks with an open, modular sovereign architecture.

---

### Slide 3: Modular System Architecture & Dataflow
- **Sequential Pipeline Flow:**
  $$\text{MPS Ingestion} \to \text{Model IR} \to \text{Validation} \to \text{Classification} \to \text{Presolve (Once)} \to \text{Engine Dispatch} \to \text{Postsolve (Once)} \to \text{Verified Solution}$$
- **Model Intermediate Representation (IR):** Explicit sparse representation of linear constraints, variable bounds, integrality types (`Continuous`, `Integer`, `Binary`), and quadratic objective matrices.
- **Strict Structural Validation:** Pre-solve validation verifies bound consistency ($lb \le ub$) and guarantees all sparse constraint indices reference valid variables.
- *Visual:* Full-width architectural flow diagram illustrating data boundaries and data-ownership handoffs between modules.

---

### Slide 4: Invertible Presolve & Dual Postsolve Restoration
- **The Single-Presolve Invariant:** Model reductions execute **exactly once** on the original model, generating an explicit audit record (`PresolveResult`). Double-presolving is strictly forbidden to prevent coordinate drift.
- **Implemented Reduction Passes:**
  - Fixed variable elimination ($lb_i = ub_i$).
  - Singleton equality row substitution ($a \cdot x_i = b$).
  - Bound tightening with provenance tracking (identifies which constraint tightened which bound).
  - Redundant constraint elimination via activity bounds.
- **Dual & Coordinate Reconstruction:** Postsolve consumes the audit metadata to restore original variable coordinates, invert bound changes, and attribute shadow prices back to original source constraints.
- **Fail-Closed Validation:** Enforces non-finite value protection (`NaN`, `±Inf`) and residual verification gates on restored solutions.
- *Visual:* Coordinate space transformation schematic showing mapping from Original Space $\to$ Reduced Space $\to$ Restored Original Space.

---

### Slide 5: Implemented Numerical Engines
- **1. PDLP Engine (`pdlp_engine/`):** First-order Primal-Dual Hybrid Gradient (PDHG) method for continuous LPs; features Ruiz diagonal equilibration, step sizes from the PDLP adaptive linesearch, and normalized gap restarts.
- **2. Dual Simplex Solver (`milp_engine/`):** Tableau-based simplex algorithm maintaining dual feasibility while driving toward primal feasibility; delivers basic feasible vertex solutions.
- **3. Branch-and-Cut Engine (`milp_engine/`):** Tree search with dynamic Gomory fractional cut generation from simplex tableau rows, most-fractional branching, and primal rounding heuristics for MILP.
- **4. Convex QP Engine (`qp_engine/`):** Alternating Direction Method of Multipliers (ADMM) with augmented KKT system factorizations for convex quadratic objectives ($f(x) = \frac{1}{2} x^T P x + q^T x + \text{offset}$).
- **5. Smooth NLP Engine (`nlp_engine/`):** Elastic SQP for `.nlp` models, reusing the QP engine. It reports first-order stationarity, not global optimality, and bypasses the affine presolve/postsolve path.
- *Visual:* Comparative matrix mapping problem families (LP, MILP, QP, NLP) to numerical engines and their underlying algorithms.

---

### Slide 6: Independent Verification Architecture
- **Verification Philosophy:** The solver is never trusted to certify its own correctness. Solutions are verified externally against the **original, unscaled model**.
- **Four-Stage Verification Pipeline:**
  1. *Dual-Reader Parse Cross-Check:* Independent C++ and Python (`mps_model.py`) parsers compare models field-by-field before solving.
  2. *Process-Level Isolation (`bench_runner`):* Each solver runs as a separate process tree with a monotonic-clock watchdog, SIGTERM-to-SIGKILL escalation, and kernel `wait4` peak RSS tracking.
  3. *Independent Solution Checker (`verify.py`):* Recomputes primal activities, bound violations, integrality, and unscaled KKT conditions against frozen tolerances ($10^{-6}$ absolute, $10^{-8}$ relative).
  4. *Separation of Status and Verdict:* Decouples solver termination claims (`optimal`, `limit_reached`) from checker proofs (`optimal_verified`, `feasible`).
- *Visual:* Mermaid flowchart showing the 4-stage pipeline from MPS file to independent verification record.

---

### Slide 7: Empirical Benchmark Evidence (Netlib LP Smoke Suite)
- **Benchmark Corpus:** A smoke suite of 8 small Netlib LP instances, not the full Netlib collection. Each solver runs in its own supervised process with a 30-second timeout and one thread. Records are in `benchmarks/results/netlib_lp.json`; the summary is `benchmarks/results/NETLIB_RESULTS.md`.
- **Empirical Results:**
  - **Dual Simplex:** Independently verifies **6 / 8** instances (`adlittle`, `afiro`, `degen2`, `recipe`, `sc50a`, `sc50b`). On `blend` and `share2b` it stops at a limit (`limit_reached`) with a feasible but clearly suboptimal point.
  - **PDLP (First-Order):** Independently verifies **6 / 8** instances (`afiro`, `blend`, `degen2`, `recipe`, `sc50a`, `share2b`). On `adlittle` and `sc50b` PDLP reports `optimal` and returns a feasible point matching the published objective, but the independent checker cannot close an optimality proof from its duals, so both are recorded as `feasible`.
  - **Together:** The two internal engines independently verify **all 8 instances** of this smoke suite (`optimal_verified`). 4 instances are verified by both engines, and each engine covers the other's 2 unverified instances. This is the union of separate single-engine runs.
  - **Reference:** The native HiGHS executable verifies 8 / 8.
- *Visual:* Smoke-suite comparison table displaying dimensions, Dual Simplex, PDLP, and native HiGHS verdicts side-by-side.

---

### Slide 8: Measured Performance & Resource Footprint
- **Measured by `bench_runner`:** Wall time comes from a monotonic clock (`std::chrono::steady_clock`) and peak RSS from `wait4(..., ru_maxrss)`. Values are means over the 8 smoke-suite instances (MB = 10⁶ bytes):
  - `optimsolver:pdlp`: mean wall time **0.090 s** | mean peak RSS **2.02 MB** (1.65–3.65 MB)
  - `optimsolver:dual_simplex`: mean wall time **0.209 s** | mean peak RSS **2.84 MB** (1.64–9.52 MB)
  - Native HiGHS reference (`highs-ds`): mean wall time **0.009 s** | mean peak RSS **4.17 MB** (3.78–5.34 MB)
- **HiGHS Methodology:** HiGHS 1.15.1 runs as its own native executable (simplex, one thread) and is measured on its own. No Python interpreter is inside the measured process.
- **Reading the Numbers:**
  - HiGHS was faster on average.
  - Mean peak RSS was lower for both internal engines on this suite, but not on every instance: Dual Simplex peaked at 9.52 MB on `degen2` against 5.34 MB for HiGHS.
  - Earlier SciPy-hosted HiGHS figures (≈ 72 MB) mostly measured the interpreter and are superseded.
- **Workload Qualification:** *Measured by the repository's isolated benchmark harness on the 8-instance Netlib smoke suite. These figures are workload-specific and are not a general performance claim.*
- *Visual:* Comparative bar charts showing peak memory (RSS) and wall times across the three configurations.

---

### Slide 9: User Experience & Dual Interface
- **Interactive Terminal Interface (`optimsolver`):** Terminal UI with ASCII mascot banner, interactive model inspection, IR statistics card, session parameter configuration, and solve summary dashboard.
- **Non-Interactive Batch Mode (`optimsolver solve <model.mps> [options]`):**
  - Engine override: `--solver <name>` (`pdlp`, `dual_simplex`, `branch_and_cut`, `qp`, `nlp` for `.nlp` models).
  - Time budgets: `--time-limit <seconds>`.
  - Solution export: `--output <file>`.
  - Machine-readable JSON export: `--json <file>`.
  - Model IR dump: `--dump-model <file>`.
  - Thread control: `--threads <n>`.
- *Visual:* Terminal screenshot collage showing the interactive main menu, model card, and batch solve output.

---

### Slide 10: Technical Integrity & Disclosed Limitations
- **Empirical Transparency:**
  - *PDLP Optimality Certification:* On 2 of the 8 smoke-suite instances (`adlittle`, `sc50b`), PDLP reports `optimal` with a feasible point that matches the published objective, but its duals do not let the independent checker close an optimality proof.
  - *Dual Simplex Limits:* On 2 of the 8 smoke-suite instances (`blend`, `share2b`), Dual Simplex stops at a limit with a feasible but suboptimal point.
  - *Benchmark Scope:* LP evidence covers an 8-instance Netlib smoke suite. The native HiGHS reference was faster on average on it.
  - *Infeasible / Unbounded Rays:* Solver engines detect infeasibility and unboundedness, but independent Farkas/improving ray certificates are not yet emitted by the checker.
  - *MILP Global Dual Bounds:* Branch-and-cut solutions are independently checked for primal feasibility and integrality (`feasible`), while global dual bound certification is actively in progress.
  - *QP Benchmark Scope:* Quadratic capabilities are verified analytically on hand-crafted KKT instances (`convex_qp.mps`, `max_qp.mps`); broad external QP benchmark corpora (e.g., Maros-Mészáros) are not claimed.
- *Visual:* High-contrast "Integrity Card" detailing current verified capabilities versus ongoing theoretical work.

---

### Slide 11: Engineering Roadmap & Sovereign Deployment
- **GPU / CUDA Acceleration (Planned):** Modular linear algebra kernels for sparse matrix-vector multiplication in PDLP and ADMM QP.
- **Infeasible Interior-Point Barrier Solver:** Interior-point methods for large continuous LPs and QPs.
- **Parallel Branch-and-Bound:** Multi-threaded tree search with concurrent node evaluation and dynamic cut pools.
- **Additional File Formats:** Native LP format (`.lp`) parser alongside existing fixed/free MPS support.
- **Sovereign Industrial Deployment:** Containerized microservice deployment for refinery blend scheduling and operational automation.
- *Visual:* Multi-phase engineering roadmap timeline from current verified foundation to GPU-accelerated enterprise deployment.

---

## Judge Takeaways

1. **Genuinely Indigenous Architecture:** Built from first principles in C++17 with clean separation across parsing, validation, invertible presolve, numerical dispatch, and postsolve.
2. **Two Internal LP Engines:** Across separate runs, Dual Simplex and PDLP together independently verify all 8 instances of an 8-instance Netlib LP smoke suite (6/8 each, 4 in common).
3. **Rigorous Independent Verification:** Backed by 64/64 automated CTest targets, dual-reader parse checks, POSIX process isolation, and external unscaled KKT verification for LP and convex QP. Infeasibility, unboundedness, and MILP dual bounds are not yet independently certified.
4. **Small Resource Footprint:** On this smoke suite the internal engines averaged about 2–2.8 MB peak RSS (maximum 9.52 MB). The native HiGHS reference averaged 4.17 MB and was faster.
5. **Scientific Honesty:** Clear demarcation between current CPU-verified solver capabilities and the future GPU acceleration roadmap.
