# Course Project: TSP Solver (Serial Baseline)

## Files

| File | Description |
|---|---|
| `TSP0.C` | Serial TSP solver (colony-based heuristic + inversion mutation) |
| `pcb442.tsp` | TSPLIB benchmark — 442 PCB drilling cities |
| `hierarchical.pdf` | Reference paper on hierarchical parallel TSP |
| `TSP0.exe` | Precompiled binary |

## Problem

Given n cities and pairwise distances, find the shortest Hamiltonian cycle visiting each city exactly once.

Benchmark `pcb442.tsp` (from TSPLIB) contains 442 cities with integer Euclidean distances (rounded).

## Algorithm

`TSP0.C` implements a population-based heuristic:

- Maintains `xColony` random permutations (tours)
- Each iteration applies **inversion mutation** to each tour: reverses a segment between two cities
- Mutation source is either random (`probab1`) or borrowed from another tour's corresponding position
- Every `xColony` iterations, parent-child tournament selection replaces worse tours
- Stops after `maxGen` generations

## Run (instructions)

Open **Developer Command Prompt for VS 2022** (or load environment via `VsDevCmd.bat`), then:

### 1. Compile

```cmd
cd course-project
cl /TC /O2 TSP0.C
```

`/TC` forces C compilation (`.C` extension otherwise treated as C++).

### 2. Run

```cmd
TSP0.exe
```

Input `pcb442.tsp` is read from current directory. Output (best distance per generation) goes to stdout:

```
init success!!!
2:720740.000000
...
2000:340177.000000
```

Detailed log is appended to `tsp0.txt`, one line per checkpoint:

```
<GenNum>  <seconds>  <distance>
```

### Tuning

Edit these variables at the top of `TSP0.C`:

| Variable | Default | Meaning |
|---|---|---|
| `xColony` | 100 | Population size |
| `probab1` | 0.02 | Random-mutation probability |
| `maxGen` | 200000 | Max generations to run |
| `NOCHANGE` | 200000 | Early-stop if no improvement |

### Known issues

- `srand` at line 119 should be `srand` but is accepted by MSVC (implicit declaration only warns)
- `%Lf` format specifier reads `long double` but `x`/`y` are `double`; works in practice on MSVC

## MPI parallelization plan (future)

Distribute `xColony` tours across processes:

1. Root scatters tour batches to workers
2. Each worker mutates + evaluates its batch independently
3. Workers send improved tours back; root performs global selection
