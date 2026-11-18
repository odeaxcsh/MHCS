# PKMC - Parallel Maximum Hyperclique Solver

A small header-only C++20 library for finding an exact maximum clique in a k-uniform hypergraph. It combines packed adjacency, core pruning, greedy coloring bounds, higher-order constraints, and a parallel worker pool with subtree donation.
The search structure is adapted from WPMC, with pruning ideas from PMC.

Vertices use zero-based, 32-bit IDs.
The order `K >= 2` is a template parameter. Each edge is a `std::array<node_t, K>` containing distinct vertices.
A clique contains every possible k-subset as an edge; for `k = 2`, this is the ordinary maximum clique problem.
Sets smaller than k are valid, so an edgeless hypergraph returns `min(n, k - 1)` vertices.

## Build

Requires CMake 3.21 or newer, a C++20 compiler, and the standard threading library.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The CMake target supplies include paths, C++20, and threading requirements;
the solver is compiled when instantiated in your project.

To use it from another CMake project:

```cmake
add_subdirectory(path/to/pkmc pkmc)
target_link_libraries(your_target PRIVATE pkmc::pkmc)
```

## Usage

```cpp
#include "clique_solver.h"

int main()
{
  const std::vector<pkmc::edge_t<3>> edges = {
    {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}, {2, 3, 4}
  };
  const pkmc::Hypergraph<3> graph(5, edges);
  const auto clique = pkmc::maximum_clique(graph, 4);
}
```

The result is one maximum clique, sorted by vertex ID: `{0, 1, 2, 3}` here.
Duplicate edges and unordered vertices within an edge are accepted. Edge size is
fixed by its type; `K < 2` is a compile-time error. Repeated vertices and
out-of-range IDs throw exceptions. `Hypergraph<2>` uses ordinary graph adjacency
without storing higher-order lookup tables.

`pkmc::CliqueSolver<K>` supports repeated sequential `search(threads)` calls and exposes `statistics()` after each search. The hypergraph must outlive the solver.
The default is four workers; zero requests one worker.
Tied maximum cliques may differ across runs.

`pkmc::SearchOptions` controls greedy seeding and reductions.
Global universal reduction is enabled by default; local universal and exchange reductions and singleton pruning are opt-in.
`reduction_budget` limits their higher-order edge lookups per search node.
