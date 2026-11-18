#pragma once

#include "hypergraph.h"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <span>
#include <vector>

namespace pkmc
{

struct SearchOptions
{
  bool global_universal_reduction = true;
  bool universal_reduction = false;
  bool exchange_reduction = false;
  bool singleton_pruning = false;
  std::size_t reduction_budget = 256;
  std::size_t greedy_starts = 16;
};

struct SearchStatistics
{
  std::uint64_t nodes = 0;
  std::uint64_t universal_forced = 0;
  std::uint64_t exchange_forced = 0;
  std::uint64_t singleton_pruned = 0;
  std::uint64_t reduction_lookups = 0;
};

std::vector<node_t> maximum_clique(
  const Hypergraph &graph, unsigned thread_count = 4, SearchOptions options = {});

class CliqueSolver
{
public:
  explicit CliqueSolver(const Hypergraph &graph, SearchOptions options = {})
    : graph_(graph), options_(options)
  {
  }

  std::vector<node_t> search(unsigned thread_count = 4);

  const SearchStatistics &statistics() const noexcept
  {
    return statistics_;
  }

private:
  struct Workspace
  {
    std::vector<node_t> cores;
    std::vector<node_t> order;
    std::vector<node_t> positions;
    std::vector<std::size_t> buckets;
    std::vector<word_t> active;
    std::vector<std::vector<node_t>> colors;
    SearchStatistics statistics;
  };

  struct Job
  {
    detail::Graph graph;
    std::vector<node_t> clique;
  };

  void update_best(std::span<const node_t> clique);

  std::vector<node_t> prepare_candidates(const detail::Graph &graph,
    std::vector<node_t> &candidates,
    std::size_t clique_size,
    Workspace &workspace);

  detail::Graph child_graph(const detail::Graph &parent,
    std::span<const node_t> candidates,
    std::span<const node_t> clique,
    node_t selected) const;

  void prune_singletons(const detail::Graph &graph,
    std::vector<node_t> &candidates,
    std::vector<node_t> &bounds,
    std::span<const node_t> clique,
    std::size_t best,
    std::size_t &budget,
    SearchStatistics &statistics) const;

  node_t try_forcing(const detail::Graph &graph,
    std::span<const node_t> candidates,
    std::span<const node_t> clique,
    std::size_t &budget,
    bool &universal,
    SearchStatistics &statistics) const;

  void branch(const detail::Graph &graph, std::vector<node_t> &clique, Workspace &workspace);

  void dispatch(detail::Graph graph, std::vector<node_t> &clique, Workspace &workspace);

  const Hypergraph &graph_;
  SearchOptions options_;
  SearchStatistics statistics_;
  bool has_global_universals_ = false;
  std::vector<node_t> cores_;
  std::vector<node_t> support_bounds_;
  std::atomic<std::size_t> best_size_{0};
  std::size_t upper_bound_ = 0;
  std::mutex best_mutex_;
  std::vector<node_t> best_clique_;

  std::mutex jobs_mutex_;
  std::condition_variable ready_;
  std::vector<Job> jobs_;
  std::atomic<unsigned> idle_{0};
  bool done_ = false;
};

}
