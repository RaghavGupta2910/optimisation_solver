# CPU vs CUDA benchmarks

* Commit: 2e87b4c-dirty (upstream main plus this branch's changes, before they were committed)
* CPU: 13th Gen Intel(R) Core(TM) i7-13645HX (20 logical processors; solver threads: all)
* GPU: NVIDIA GeForce RTX 5050 Laptop GPU, 12.0, 610.74, 8151 MiB
* CUDA Toolkit: 13.4
* OS: Microsoft Windows NT 10.0.26200.0
* Method: Release build, double precision, 1 untimed warm-up + 5 timed end-to-end solves per backend; median reported; CPU and CUDA run back to back per instance
* Speedup = CPU median / CUDA median; objective error = abs(CPU objective - CUDA objective)

## PDLP: CPU vs CUDA

Generated feasible LPs (pdlp_bench, seed 12345): rows x 2*rows columns, 10 nonzeros per row, iteration limit 200000, default tolerances (1e-6), polishing off.

| Instance | NNZ | CPU (s) | CUDA (s) | Speedup | Objective Error | CPU / CUDA iterations | CPU / CUDA status | CUDA setup (s) | CPU min-max (s) | CUDA min-max (s) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: |
| 1000x2000x10 | 9,977 | 0.0806 | 0.454 | 0.18x | 1.61e-05 | 6000 / 6000 | optimal / optimal | 0.0006 | 0.0794-0.0906 | 0.451-0.484 |
| 10000x20000x10 | 99,976 | 0.793 | 2.251 | 0.35x | 1.57e-05 | 11200 / 11200 | optimal / optimal | 0.0042 | 0.645-7.323 | 1.579-2.452 |
| 25000x50000x10 | 249,972 | 8.080 | 2.390 | 3.38x | 8.51e-06 | 7800 / 7800 | optimal / optimal | 0.0065 | 1.043-8.343 | 1.684-2.479 |
| 50000x100000x10 | 499,974 | 11.274 | 2.781 | 4.05x | 3.56e-05 | 8000 / 8000 | optimal / optimal | 0.0063 | 2.224-16.370 | 2.677-3.082 |
| 75000x150000x10 | 749,971 | 5.813 | 5.929 | 0.98x | 2.21e-04 | 9000 / 9000 | optimal / optimal | 0.0149 | 5.148-29.732 | 5.908-5.951 |
| 100000x200000x10 | 999,977 | 40.438 | 7.657 | 5.28x | 1.26e-05 | 9000 / 9000 | optimal / optimal | 0.0159 | 39.009-40.956 | 7.133-7.702 |
| 150000x300000x10 | 1,499,978 | 59.545 | 10.667 | 5.58x | 2.72e-04 | 8600 / 8800 | optimal / optimal | 0.0315 | 51.920-61.052 | 10.498-10.695 |
| 200000x400000x10 | 1,999,974 | 80.831 | 17.755 | 4.55x | 5.64e-05 | 8600 / 8600 | optimal / optimal | 0.0208 | 63.262-82.612 | 10.657-29.563 |
| 500000x1000000x10 | 4,999,978 | 149.773 | 37.739 | 3.97x | 7.58e-05 | 8800 / 8800 | optimal / optimal | 0.0886 | 84.948-352.515 | 37.596-37.927 |

## QP: CPU vs hybrid CUDA

Generated banded convex QPs (`qp_bench sweep`, seed 2024): n variables, n constraint rows with 8 nonzeros each, tridiagonal P, default ADMM options. The CUDA backend is hybrid: the KKT factorisation and solves stay on the CPU.

| Instance | NNZ (P+A) | CPU (s) | Hybrid CUDA (s) | Speedup | CUDA KKT solve (s) | CUDA setup (s) | H2D / D2H (MB) | CPU / CUDA iterations | CPU / CUDA status | Objective Error |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: |
| 2000x8 | 18,869 | 0.0230 | 0.0514 | 0.45x | 0.0066 | 0.0012 | 2.7 / 2.5 | 150 / 150 | Optimal / Optimal | 0.00e+00 |
| 5000x8 | 47,209 | 0.125 | 0.152 | 0.82x | 0.0196 | 0.0014 | 6.8 / 6.2 | 150 / 150 | Optimal / Optimal | 4.00e-14 |
| 10000x8 | 94,375 | 0.448 | 0.500 | 0.90x | 0.0423 | 0.0020 | 13.6 / 12.4 | 150 / 150 | Optimal / Optimal | 5.70e-13 |
| 20000x8 | 189,035 | 5.514 | 5.771 | 0.96x | 0.749 | 0.0079 | 42.5 / 41.2 | 250 / 250 | Optimal / Optimal | 8.00e-13 |
| 30000x8 | 283,425 | 11.816 | 11.579 | 1.02x | 1.165 | 0.0099 | 63.7 / 61.8 | 250 / 250 | Optimal / Optimal | 1.03e-12 |
