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
- *Visual:* MRPL & SIH partner branding, team grid, and repository badges (C++17, CMake 3.20+, CTest 60/60, MIT).

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
- **1. PDLP Engine (`pdlp_engine/`):** First-order Primal-Dual Hybrid Gradient (PDHG) method for continuous LPs; features Ruiz diagonal equilibration, adaptive Malitsky-Pock step sizes, and normalized gap restarts.
- **2. Dual Simplex Solver (`milp_engine/`):** Tableau-based simplex algorithm maintaining dual feasibility while driving toward primal feasibility; delivers basic feasible vertex solutions.
- **3. Branch-and-Cut Engine (`milp_engine/`):** Tree search with dynamic Gomory fractional cut generation from simplex tableau rows, most-fractional branching, and primal rounding heuristics for MILP.
- **4. Convex QP Engine (`qp_engine/`):** Alternating Direction Method of Multipliers (ADMM) with augmented KKT system factorizations for convex quadratic objectives ($f(x) = \frac{1}{2} x^T P x + q^T x + \text{offset}$).
- *Visual:* Comparative matrix mapping problem families (LP, MILP, QP) to numerical engines and their underlying algorithms.

---

### Slide 6: Independent Verification Architecture
- **Verification Philosophy:** The solver is never trusted to certify its own correctness. Solutions are verified externally against the **original, unscaled model**.
- **Four-Stage Verification Pipeline:**
  1. *Dual-Reader Parse Cross-Check:* Independent C++ and Python (`mps_model.py`) parsers compare models field-by-field before solving.
  2. *Process-Level Isolation (`bench_runner`):* Separate process groups with POSIX monotonic watchdogs, SIGTERM-to-SIGKILL escalation, and kernel `wait4` peak RSS tracking.
  3. *Independent Solution Checker (`verify.py`):* Recomputes primal activities, bound violations, integrality, and unscaled KKT conditions against frozen tolerances ($10^{-6}$ absolute, $10^{-8}$ relative).
  4. *Separation of Status and Verdict:* Decouples solver termination claims (`optimal`, `limit_reached`) from checker proofs (`optimal_verified`, `feasible`).
- *Visual:* Mermaid flowchart showing the 4-stage pipeline from MPS file to independent verification record.

---

### Slide 7: Empirical Benchmark Evidence (Netlib LP Suite)
- **Benchmark Corpus:** Standard Netlib LP smoke benchmark suite (8 instances) evaluated with a 60-second time budget under isolated process supervision (`benchmarks/results/netlib_smoke.json`).
- **Empirical Results:**
  - **Dual Simplex:** Independently verifies **6 / 8** instances (`adlittle`, `afiro`, `degen2`, `recipe`, `sc50a`, `sc50b`); encounters iteration limits on cycling/dense models (`blend`, `share2b`).
  - **PDLP (First-Order):** Independently verifies **6 / 8** instances (`afiro`, `blend`, `degen2`, `recipe`, `sc50a`, `share2b`); encounters numerical stalls at the tight relative tolerance boundary on `adlittle` and `sc50b`.
  - **Combined Portfolio Coverage:** **8 / 8 instances (100%) independently verified** (`optimal_verified`). Dual Simplex solves the instances PDLP misses, and PDLP solves the instances Dual Simplex stalls on.
- *Visual:* Netlib smoke comparison table displaying dimensions, Dual Simplex, PDLP, and HiGHS verdicts side-by-side.

---

### Slide 8: Measured Performance & Resource Footprint
- **Kernel-Measured Resource Footprint:** Monitored via POSIX `CLOCK_MONOTONIC` and `wait4(..., ru_maxrss)`:
  - `optimsolver:pdlp`: Mean wall time **0.056 s** | Mean peak RSS **1.98 MB**
  - `optimsolver:dual_simplex`: Mean wall time **0.146 s** | Mean peak RSS **2.83 MB**
  - HiGHS Reference Process (`highs-ds`): Mean wall time **0.340 s** | Mean peak RSS **71.82 MB**
- **Low Memory Overhead:** Native C++ implementation operates in ~1.6–2.8 MB RSS, exhibiting low runtime overhead compared to external wrapper processes.
- **Workload Qualification:** *Measured by the repository's isolated benchmark harness on the current 8-instance Netlib smoke set. These figures are workload-specific and are not a general performance claim.*
- *Visual:* Comparative bar charts showing peak memory (RSS) and solve times across engines.

---

### Slide 9: User Experience & Dual Interface
- **Interactive Terminal Interface (`optimsolver`):** Terminal UI with ASCII mascot banner, interactive model inspection, IR statistics card, session parameter configuration, and solve summary dashboard.
- **Non-Interactive Batch Mode (`optimsolver solve <model.mps> [options]`):**
  - Engine override: `--solver <name>` (`pdlp`, `dual_simplex`, `branch_and_cut`, `qp`).
  - Time budgets: `--time-limit <seconds>`.
  - Solution export: `--output <file>`.
  - Machine-readable JSON export: `--json <file>`.
  - Model IR dump: `--dump-model <file>`.
  - Thread control: `--threads <n>`.
- *Visual:* Terminal screenshot collage showing the interactive main menu, model card, and batch solve output.

---

### Slide 10: Technical Integrity & Disclosed Limitations
- **Empirical Transparency:**
  - *PDLP Boundary Tuning:* PDLP encounters convergence stalls at tight relative tolerances ($10^{-8}$) on ill-conditioned bases (2/8 Netlib cases).
  - *Infeasible / Unbounded Rays:* Solver engines detect infeasibility and unboundedness, but independent Farkas/improving ray certificates are not yet emitted by the checker.
  - *MILP Global Dual Bounds:* Branch-and-cut guarantees integer feasibility and cost evaluation (`feasible`), while global dual bound certification is actively in progress.
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
2. **Complementary Engine Portfolio:** Combined internal engines (Dual Simplex + PDLP) independently verify 100% (8/8) of the standard Netlib LP smoke benchmark suite.
3. **Rigorous Independent Verification:** Backed by 60/60 automated CTest targets, dual-reader parse checks, POSIX process isolation, and external unscaled KKT verification.
4. **Lightweight Resource Footprint:** Native C++ implementation achieves ~2 MB peak RSS on Netlib instances, providing an efficient foundation for industrial edge and server deployments.
5. **Scientific Honesty:** Clear demarcation between current CPU-verified solver capabilities and the future GPU acceleration roadmap.
