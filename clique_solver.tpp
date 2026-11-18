#pragma once

#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <thread>
#include <utility>

namespace pkmc
{
namespace detail
{

inline constexpr node_t None = std::numeric_limits<node_t>::max();

inline bool binomial_exceeds(std::size_t n, std::size_t k, std::size_t limit)
{
  if(k > n) {
    return false;
  }
  k = std::min(k, n - k);
  std::size_t value = 1;
  for(std::size_t i = 1; i <= k; ++i) {
    std::size_t numerator = n - k + i;
    std::size_t denominator = i;
    const auto divisor = std::gcd(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;
    value /= denominator; // The reduced denominator divides the previous binomial coefficient.
    if(value > limit / numerator) {
      return true;
    }
    value *= numerator;
  }
  return value > limit;
}

}

template <std::size_t K>
void CliqueSolver<K>::update_best(std::span<const node_t> clique)
{
  if(clique.size() <= best_size_.load(std::memory_order_relaxed)) {
    return;
  }
  std::lock_guard lock(best_mutex_);
  if(clique.size() <= best_size_.load(std::memory_order_relaxed)) {
    return;
  }
  best_clique_.assign(clique.begin(), clique.end());
  std::sort(best_clique_.begin(), best_clique_.end());
  best_size_.store(clique.size(), std::memory_order_relaxed);
}

template <std::size_t K>
std::vector<node_t> CliqueSolver<K>::search(unsigned thread_count)
{
  statistics_ = {};
  best_size_ = 0;
  best_clique_.clear();
  jobs_.clear();
  idle_ = 0;
  done_ = false;
  const auto &projection = graph_.projection();

  std::vector<node_t> baseline(std::min(graph_.size(), K - 1));
  std::iota(baseline.begin(), baseline.end(), node_t{0});
  update_best(baseline);
  if(!graph_.edges().empty()) {
    update_best(graph_.edges().front());
  }
  if(graph_.edges().empty() || best_size_ == graph_.size()) {
    return best_clique_;
  }

  support_bounds_.resize(graph_.size());
  for(std::size_t v = 0; v < graph_.size(); ++v) {
    if constexpr(K == 2) {
      support_bounds_[v] = graph_.incidence(v) + 1;
    } else {
      std::size_t low = K - 1;
      std::size_t high = graph_.size();
      while(low < high) {
        const auto middle = low + (high - low + 1) / 2;
        if(detail::binomial_exceeds(middle - 1, K - 1, graph_.incidence(v))) {
          high = middle - 1;
        } else {
          low = middle;
        }
      }
      support_bounds_[v] = low;
    }
  }
  has_global_universals_ =
    std::ranges::find(support_bounds_, graph_.size()) != support_bounds_.end();

  Workspace initial;
  std::vector<node_t> candidates(graph_.size());
  std::iota(candidates.begin(), candidates.end(), node_t{0});
  auto bounds = prepare_candidates(projection, candidates, 0, initial);
  cores_ = initial.cores;
  if(candidates.empty()) {
    return best_clique_;
  }
  auto order = initial.order;
  upper_bound_ =
    std::min<std::size_t>(*std::max_element(cores_.begin(), cores_.end()) + 1, bounds.back());

  for(std::size_t start = 0; start < std::min(options_.greedy_starts, order.size()); ++start) {
    std::vector<node_t> greedy{projection.vertex(order[start])};
    for(node_t local : order) {
      const auto vertex = projection.vertex(local);
      if(vertex != greedy.front() &&
         std::ranges::all_of(greedy, [&](node_t v) { return projection.is_edge(vertex, v); })) {
        if constexpr(K == 2) {
          greedy.push_back(vertex);
        } else if(graph_.can_extend(greedy, vertex)) {
          greedy.push_back(vertex);
        }
      }
    }
    update_best(greedy);
    if(best_size_ >= upper_bound_) {
      return best_clique_;
    }
  }

  std::vector<node_t> global_color(graph_.size());
  for(std::size_t i = 0; i < candidates.size(); ++i) {
    global_color[candidates[i]] = bounds[i];
  }
  std::vector<std::size_t> rank(graph_.size(), order.size());
  for(std::size_t i = 0; i < order.size(); ++i) {
    rank[order[i]] = i;
  }

  thread_count = std::max<std::size_t>(1, std::min<std::size_t>(thread_count, order.size()));
  std::vector<SearchStatistics> worker_statistics(thread_count);
  std::atomic<std::size_t> next{0};
  auto work = [&](unsigned id) {
    Workspace workspace = id == 0 ? std::move(initial) : Workspace{};
    std::vector<node_t> clique;
    std::vector<node_t> remaining;
    std::vector<std::size_t> color_seen(bounds.back() + 1, order.size());
    auto root = [&](std::size_t index) {
      const auto best = best_size_.load(std::memory_order_relaxed);
      const auto selected = order[index];
      if(best >= upper_bound_ || cores_[selected] < best || support_bounds_[selected] <= best) {
        return;
      }
      remaining.clear();
      node_t distinct_colors = 0;
      for(node_t vertex : order) {
        if(rank[vertex] > index && cores_[vertex] >= best && support_bounds_[vertex] > best &&
           projection.is_edge(selected, vertex)) {
          remaining.push_back(vertex);
          const auto color = global_color[vertex];
          if(color_seen[color] != index) {
            color_seen[color] = index;
            ++distinct_colors;
          }
        }
      }
      if(remaining.size() + 1 <= best || distinct_colors + 1 <= best) {
        return;
      }
      auto child = child_graph(projection, remaining, {}, selected);
      clique.assign(1, projection.vertex(selected));
      branch(child, clique, workspace);
    };
    for(;;) {
      const auto index = next.fetch_add(1, std::memory_order_relaxed);
      if(index < order.size()) {
        root(index);
        continue;
      }
      std::unique_lock lock(jobs_mutex_);
      ++idle_;
      if(idle_ == thread_count && jobs_.empty()) {
        done_ = true;
        ready_.notify_all();
      }
      ready_.wait(lock, [&] { return done_ || !jobs_.empty(); });
      --idle_;
      if(done_) {
        break;
      }
      auto job = std::move(jobs_.back());
      jobs_.pop_back();
      lock.unlock();
      branch(job.graph, job.clique, workspace);
    }
    worker_statistics[id] = workspace.statistics;
  };
  {
    std::vector<std::jthread> workers;
    workers.reserve(thread_count - 1);
    for(unsigned id = 1; id < thread_count; ++id) {
      workers.emplace_back(work, id);
    }
    work(0);
  }
  for(const auto &s : worker_statistics) {
    statistics_.nodes += s.nodes;
    statistics_.universal_forced += s.universal_forced;
    statistics_.exchange_forced += s.exchange_forced;
    statistics_.singleton_pruned += s.singleton_pruned;
    statistics_.reduction_lookups += s.reduction_lookups;
  }
  std::sort(best_clique_.begin(), best_clique_.end());
  return best_clique_;
}

template <std::size_t K>
std::vector<node_t> CliqueSolver<K>::prepare_candidates(const detail::Graph &graph,
  std::vector<node_t> &candidates,
  std::size_t clique_size,
  Workspace &workspace)
{
  const auto best = best_size_.load(std::memory_order_relaxed);
  const auto minimum_size = best >= clique_size ? best - clique_size + 1 : 1;
  if(candidates.size() < minimum_size) {
    candidates.clear();
    return {};
  }
  auto &active = workspace.active;
  active.assign((graph.size() + WordBits - 1) / WordBits, 0);
  for(node_t v : candidates) {
    active[v / WordBits] |= word_t{1} << (v % WordBits);
  }
  auto &cores = workspace.cores;
  cores.assign(graph.size(), 0);
  node_t maximum_degree = 0;
  for(node_t v : candidates) {
    const auto row = graph.adjacency(v);
    for(std::size_t w = 0; w < active.size(); ++w) {
      cores[v] += std::popcount(row[w] & active[w]);
    }
    maximum_degree = std::max(maximum_degree, cores[v]);
  }

  auto &buckets = workspace.buckets;
  buckets.assign(std::size_t(maximum_degree) + 1, 0);
  for(node_t v : candidates) {
    ++buckets[cores[v]];
  }
  std::size_t start = 0;
  for(auto &bucket : buckets) {
    const auto count = bucket;
    bucket = start;
    start += count;
  }
  auto &positions = workspace.positions;
  positions.resize(graph.size());
  auto &order = workspace.order;
  order.resize(candidates.size());
  for(node_t v : candidates) {
    positions[v] = buckets[cores[v]]++;
    order[positions[v]] = v;
  }
  for(std::size_t d = buckets.size() - 1; d > 0; --d) {
    buckets[d] = buckets[d - 1];
  }
  buckets[0] = 0;
  for(std::size_t i = 0; i < order.size(); ++i) {
    const auto v = order[i];
    active[v / WordBits] &= ~(word_t{1} << (v % WordBits));
    const auto row = graph.adjacency(v);
    for(std::size_t w = 0; w < active.size(); ++w) {
      auto neighbors = row[w] & active[w];
      while(neighbors) {
        const node_t u = w * WordBits + std::countr_zero(neighbors);
        neighbors &= neighbors - 1;
        if(cores[u] <= cores[v]) {
          continue;
        }
        const auto degree = cores[u];
        const auto position = positions[u];
        const auto first_position = buckets[degree];
        const auto first = order[first_position];
        std::swap(order[position], order[first_position]);
        positions[u] = first_position;
        positions[first] = position;
        ++buckets[degree];
        --cores[u];
      }
    }
  }
  std::reverse(order.begin(), order.end());
  std::erase_if(order, [&](node_t v) {
    return cores[v] + 1 < minimum_size || support_bounds_[graph.vertex(v)] <= best;
  });
  if(order.size() < minimum_size) {
    candidates.clear();
    return {};
  }

  auto &colors = workspace.colors;
  colors.resize(order.size());
  std::size_t used_colors = 0;
  for(node_t v : order) {
    std::size_t color = 0;
    while(color < used_colors &&
          std::ranges::any_of(colors[color], [&](node_t u) { return graph.is_edge(v, u); })) {
      ++color;
    }
    if(color == used_colors) {
      colors[used_colors++].clear();
    }
    colors[color].push_back(v);
  }
  candidates.clear();
  std::vector<node_t> bounds;
  bounds.reserve(order.size());
  for(std::size_t color = 0; color < used_colors; ++color) {
    candidates.insert(candidates.end(), colors[color].begin(), colors[color].end());
    bounds.insert(bounds.end(), colors[color].size(), color + 1);
  }
  return bounds;
}

template <std::size_t K>
detail::Graph CliqueSolver<K>::child_graph(const detail::Graph &parent,
  std::span<const node_t> candidates,
  std::span<const node_t> clique,
  node_t selected) const
{
  std::vector<node_t> vertices;
  vertices.reserve(candidates.size());
  for(node_t local : candidates) {
    vertices.push_back(parent.vertex(local));
  }
  detail::Graph child(std::move(vertices));
  for(std::size_t i = 0; i < candidates.size(); ++i) {
    for(std::size_t j = i + 1; j < candidates.size(); ++j) {
      if(parent.is_edge(candidates[i], candidates[j])) {
        child.add_edge(i, j);
      }
    }
  }
  if constexpr(K > 2) {
    if(child.size() < 2) {
      return child;
    }
    edge_t<K - 3> base{};
    detail::combinations(clique, base, [&] {
      edge_t<K - 1> face{};
      std::copy(base.begin(), base.end(), face.begin());
      face[K - 3] = parent.vertex(selected);
      for(std::size_t i = 0; i < child.size(); ++i) {
        face.back() = child.vertex(i);
        auto sorted_face = face;
        std::sort(sorted_face.begin(), sorted_face.end());
        const auto *neighbors = graph_.link(sorted_face);
        auto row = child.adjacency(i);
        if(!neighbors) {
          std::fill(row.begin(), row.end(), 0);
          continue;
        }
        for(std::size_t w = 0; w < row.size(); ++w) {
          auto bits = row[w];
          while(bits) {
            const auto offset = std::countr_zero(bits);
            const node_t j = w * WordBits + offset;
            if(!neighbors->contains(child.vertex(j))) {
              row[w] &= ~(word_t{1} << offset);
            }
            bits &= bits - 1;
          }
        }
      }
      return true;
    });
  }
  return child;
}

template <std::size_t K>
void CliqueSolver<K>::prune_singletons(const detail::Graph &graph,
  std::vector<node_t> &candidates,
  std::vector<node_t> &bounds,
  std::span<const node_t> clique,
  std::size_t best,
  std::size_t &budget,
  SearchStatistics &statistics) const
{
  std::vector<node_t> singletons;
  for(std::size_t i = 0; i < candidates.size(); ++i) {
    if((i == 0 || bounds[i - 1] != bounds[i]) &&
       (i + 1 == bounds.size() || bounds[i + 1] != bounds[i])) {
      singletons.push_back(candidates[i]);
    }
  }
  const auto slack = clique.size() + bounds.back() - best - 1;
  if(singletons.size() <= slack) {
    return;
  }
  std::vector<word_t> singleton_mask((graph.size() + WordBits - 1) / WordBits, 0);
  for(node_t v : singletons) {
    singleton_mask[v / WordBits] |= word_t{1} << (v % WordBits);
  }
  std::vector<bool> removed(graph.size(), false);
  std::vector<bool> used;
  if constexpr(K > 2) {
    used.resize(graph.size());
  }
  for(node_t x : candidates) {
    std::size_t packed = 0;
    const auto row = graph.adjacency(x);
    for(std::size_t w = 0; w < singleton_mask.size(); ++w) {
      auto conflicts = singleton_mask[w] & ~row[w];
      if(w == x / WordBits) {
        conflicts &= ~(word_t{1} << (x % WordBits));
      }
      packed += std::popcount(conflicts);
    }
    if constexpr(K > 2) {
      if(packed <= slack && budget) {
        std::fill(used.begin(), used.end(), false);
        for(node_t v : singletons) {
          used[v] = v == x || !graph.is_edge(x, v);
        }
      }
      for(std::size_t r = 2; packed <= slack && budget && r < K && r <= singletons.size(); ++r) {
        if(K - 1 - r > clique.size()) {
          continue;
        }
        bool insufficient = false;
        while(packed <= slack && budget) {
          std::vector<node_t> available;
          for(node_t v : singletons) {
            if(!used[v]) {
              available.push_back(v);
            }
          }
          if(available.size() / r <= slack - packed) {
            insufficient = true;
            break;
          }
          bool found = false;
          edge_t<K - 1> subset{};
          detail::combinations(available, std::span(subset).first(r), [&] {
            bool missing = false;
            edge_t<K - 1> face{};
            detail::combinations(clique, std::span(face).first(K - 1 - r), [&] {
              if(!budget) {
                return false;
              }
              --budget;
              ++statistics.reduction_lookups;
              auto sorted_face = face;
              for(std::size_t i = 0; i < r; ++i) {
                sorted_face[K - 1 - r + i] = graph.vertex(subset[i]);
              }
              std::sort(sorted_face.begin(), sorted_face.end());
              const auto *neighbors = graph_.link(sorted_face);
              missing = !neighbors || !neighbors->contains(graph.vertex(x));
              return !missing;
            });
            if(missing) {
              for(node_t v : std::span(subset).first(r)) {
                used[v] = true;
              }
              ++packed;
              found = true;
            }
            return !found && budget;
          });
          if(!found) {
            break;
          }
        }
        if(insufficient) {
          break;
        }
      }
    }
    if(packed > slack) {
      removed[x] = true;
      ++statistics.singleton_pruned;
    }
  }
  std::size_t retained = 0;
  node_t color = 0;
  node_t previous = 0;
  for(std::size_t i = 0; i < candidates.size(); ++i) {
    if(removed[candidates[i]]) {
      continue;
    }
    if(bounds[i] != previous) {
      previous = bounds[i];
      ++color;
    }
    candidates[retained] = candidates[i];
    bounds[retained++] = color;
  }
  candidates.resize(retained);
  bounds.resize(retained);
}

template <std::size_t K>
node_t CliqueSolver<K>::try_forcing(const detail::Graph &graph,
  std::span<const node_t> candidates,
  std::span<const node_t> clique,
  std::size_t &budget,
  bool &universal,
  SearchStatistics &statistics) const
{
  if(has_global_universals_ &&
     (options_.global_universal_reduction || options_.universal_reduction)) {
    for(node_t v : candidates) {
      if(support_bounds_[graph.vertex(v)] == graph_.size()) {
        universal = true;
        ++statistics.universal_forced;
        return v;
      }
    }
  }
  if(!options_.universal_reduction && !options_.exchange_reduction) {
    return detail::None;
  }
  std::vector<word_t> active((graph.size() + WordBits - 1) / WordBits, 0);
  for(node_t v : candidates) {
    active[v / WordBits] |= word_t{1} << (v % WordBits);
  }
  for(node_t v : candidates) {
    node_t exception = detail::None;
    bool multiple = false;
    const auto row = graph.adjacency(v);
    for(std::size_t w = 0; w < active.size(); ++w) {
      auto conflicts = active[w] & ~row[w];
      if(w == v / WordBits) {
        conflicts &= ~(word_t{1} << (v % WordBits));
      }
      if(!conflicts) {
        continue;
      }
      if(exception != detail::None || !options_.exchange_reduction ||
         std::popcount(conflicts) > 1) {
        multiple = true;
        break;
      }
      exception = w * WordBits + std::countr_zero(conflicts);
    }
    if(multiple) {
      continue;
    }
    if constexpr(K == 2) {
      universal = exception == detail::None;
      if(universal && !options_.universal_reduction) {
        continue;
      }
    } else {
      if(!budget) {
        continue;
      }
      const auto excluded = options_.exchange_reduction ? 2 : 1;
      if(detail::binomial_exceeds(clique.size() + candidates.size() - excluded, K - 1,
           graph_.incidence(graph.vertex(v)))) {
        continue;
      }
      std::vector<node_t> others;
      for(node_t u : candidates) {
        if(u != v && u != exception) {
          others.push_back(u);
        }
      }
      edge_t<K - 1> common{};
      std::size_t common_size = 0;
      bool missing = false;
      bool complete = true;
      for(std::size_t r = 2; complete && r < K && r <= others.size(); ++r) {
        if(K - 1 - r > clique.size()) {
          continue;
        }
        edge_t<K - 1> subset{};
        complete = detail::combinations(others, std::span(subset).first(r), [&] {
          edge_t<K - 1> face{};
          return detail::combinations(clique, std::span(face).first(K - 1 - r), [&] {
            if(!budget) {
              return false;
            }
            --budget;
            ++statistics.reduction_lookups;
            auto sorted_face = face;
            for(std::size_t i = 0; i < r; ++i) {
              sorted_face[K - 1 - r + i] = graph.vertex(subset[i]);
            }
            std::sort(sorted_face.begin(), sorted_face.end());
            const auto *neighbors = graph_.link(sorted_face);
            if(neighbors && neighbors->contains(graph.vertex(v))) {
              return true;
            }
            if(exception != detail::None || !options_.exchange_reduction) {
              return false;
            }
            if(!missing) {
              std::copy_n(subset.begin(), r, common.begin());
              common_size = r;
              missing = true;
            } else {
              const auto selected = std::span(subset).first(r);
              common_size =
                std::remove_if(common.begin(), common.begin() + common_size,
                  [&](node_t u) { return std::ranges::find(selected, u) == selected.end(); }) -
                common.begin();
            }
            return common_size != 0;
          });
        });
      }
      if(!complete) {
        continue;
      }
      universal = exception == detail::None && !missing;
      if(universal && !options_.universal_reduction) {
        continue;
      }
    }
    if(universal) {
      ++statistics.universal_forced;
    } else {
      ++statistics.exchange_forced;
    }
    return v;
  }
  return detail::None;
}

template <std::size_t K>
void CliqueSolver<K>::dispatch(
  detail::Graph graph, std::vector<node_t> &clique, Workspace &workspace)
{
  if(idle_.load(std::memory_order_relaxed)) {
    std::unique_lock lock(jobs_mutex_);
    if(jobs_.size() < idle_.load(std::memory_order_relaxed)) {
      jobs_.push_back({std::move(graph), clique});
      lock.unlock();
      ready_.notify_one();
      return;
    }
  }
  branch(graph, clique, workspace);
}

template <std::size_t K>
void CliqueSolver<K>::branch(
  const detail::Graph &graph, std::vector<node_t> &clique, Workspace &workspace)
{
  ++workspace.statistics.nodes;
  update_best(clique);
  auto best = best_size_.load(std::memory_order_relaxed);
  if(best >= upper_bound_ || clique.size() + graph.size() <= best) {
    return;
  }
  std::vector<node_t> candidates(graph.size());
  std::iota(candidates.begin(), candidates.end(), node_t{0});
  auto bounds = prepare_candidates(graph, candidates, clique.size(), workspace);
  const auto entry_size = clique.size();
  auto reduction_budget = options_.reduction_budget;
  while(!candidates.empty()) {
    best = best_size_.load(std::memory_order_relaxed);
    if(best >= upper_bound_ || clique.size() + bounds.back() <= best) {
      break;
    }
    if(options_.singleton_pruning) {
      prune_singletons(
        graph, candidates, bounds, clique, best, reduction_budget, workspace.statistics);
      if(candidates.empty() || clique.size() + bounds.back() <= best) {
        break;
      }
    }
    bool universal = false;
    const auto forced =
      options_.universal_reduction || options_.exchange_reduction ||
          (options_.global_universal_reduction && has_global_universals_)
        ? try_forcing(graph, candidates, clique, reduction_budget, universal, workspace.statistics)
        : detail::None;
    if(forced != detail::None) {
      if(universal) {
        const auto index = std::ranges::find(candidates, forced) - candidates.begin();
        const auto color = bounds[index];
        candidates.erase(candidates.begin() + index);
        bounds.erase(bounds.begin() + index);
        for(auto &bound : bounds) {
          if(bound > color) {
            --bound;
          }
        }
        clique.push_back(graph.vertex(forced));
        update_best(clique);
        continue;
      }
      std::vector<node_t> remaining;
      for(node_t v : candidates) {
        if(v != forced && graph.is_edge(forced, v)) {
          remaining.push_back(v);
        }
      }
      auto child = child_graph(graph, remaining, clique, forced);
      clique.push_back(graph.vertex(forced));
      dispatch(std::move(child), clique, workspace);
      clique.pop_back();
      break;
    }
    const auto selected = candidates.back();
    candidates.pop_back();
    bounds.pop_back();
    std::vector<node_t> remaining;
    node_t distinct_colors = 0;
    node_t previous_color = 0;
    for(std::size_t i = 0; i < candidates.size(); ++i) {
      const auto v = candidates[i];
      if(graph.is_edge(selected, v) && cores_[graph.vertex(v)] >= best &&
         support_bounds_[graph.vertex(v)] > best) {
        remaining.push_back(v);
        if(bounds[i] != previous_color) {
          previous_color = bounds[i];
          ++distinct_colors;
        }
      }
    }
    if(remaining.empty()) {
      clique.push_back(graph.vertex(selected));
      update_best(clique);
      clique.pop_back();
    } else if(clique.size() + 1 + distinct_colors > best) {
      auto child = child_graph(graph, remaining, clique, selected);
      clique.push_back(graph.vertex(selected));
      dispatch(std::move(child), clique, workspace);
      clique.pop_back();
    }
  }
  clique.resize(entry_size);
}

template <std::size_t K>
std::vector<node_t> maximum_clique(
  const Hypergraph<K> &graph, unsigned thread_count, SearchOptions options)
{
  return CliqueSolver<K>(graph, options).search(thread_count);
}

}
