# CPU vs CUDA benchmarks

* Commit: 2e87b4c-dirty (upstream main plus this branch's changes, before they were committed)
* CPU: 13th Gen Intel(R) Core(TM) i7-13645HX (20 logical processors; solver threads: 12)
* GPU: NVIDIA GeForce RTX 5050 Laptop GPU, 12.0, 610.74, 8151 MiB
* CUDA Toolkit: 13.4
* OS: Microsoft Windows NT 10.0.26200.0
* Method: Release build, double precision, 1 untimed warm-up + 5 timed end-to-end solves per backend; median reported; CPU and CUDA run back to back per instance
* Speedup = CPU median / CUDA median; objective error = abs(CPU objective - CUDA objective)

## PDLP: CPU vs CUDA

Generated feasible LPs (pdlp_bench, seed 12345): rows x 2*rows columns, 10 nonzeros per row, iteration limit 200000, default tolerances (1e-6), polishing off.

| Instance | NNZ | CPU (s) | CUDA (s) | Speedup | Objective Error | CPU / CUDA iterations | CPU / CUDA status | CUDA setup (s) | CPU min-max (s) | CUDA min-max (s) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: |
| 1000x2000x10 | 9,977 | 0.0742 | 0.359 | 0.21x | 1.61e-05 | 6000 / 6000 | optimal / optimal | 0.0006 | 0.0709-0.0775 | 0.351-0.368 |
| 10000x20000x10 | 99,976 | 0.464 | 1.262 | 0.37x | 2.15e-05 | 11200 / 11200 | optimal / optimal | 0.0018 | 0.461-0.471 | 1.254-1.294 |
| 25000x50000x10 | 249,972 | 0.975 | 1.493 | 0.65x | 6.92e-05 | 7800 / 7800 | optimal / optimal | 0.0023 | 0.830-1.065 | 1.464-1.530 |
| 50000x100000x10 | 499,974 | 2.313 | 2.529 | 0.91x | 5.48e-05 | 8000 / 8000 | optimal / optimal | 0.0044 | 2.288-2.331 | 2.499-2.539 |
| 75000x150000x10 | 749,971 | 4.764 | 4.073 | 1.17x | 2.62e-05 | 9000 / 9000 | optimal / optimal | 0.0062 | 4.711-4.838 | 4.021-4.183 |
| 100000x200000x10 | 999,977 | 8.178 | 5.480 | 1.49x | 4.29e-06 | 9000 / 9000 | optimal / optimal | 0.0084 | 8.141-8.497 | 5.389-5.485 |
| 150000x300000x10 | 1,499,978 | 16.746 | 7.667 | 2.18x | 1.49e-04 | 8600 / 8800 | optimal / optimal | 0.0184 | 16.058-18.544 | 7.611-7.767 |
| 200000x400000x10 | 1,999,974 | 24.762 | 9.845 | 2.52x | 2.62e-05 | 8600 / 8600 | optimal / optimal | 0.0137 | 24.629-24.906 | 9.731-10.005 |
| 500000x1000000x10 | 4,999,978 | 75.848 | 26.256 | 2.89x | 1.66e-04 | 8800 / 8800 | optimal / optimal | 0.0348 | 75.078-79.043 | 25.946-32.023 |
