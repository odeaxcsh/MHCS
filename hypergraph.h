#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pkmc
{

using node_t = std::uint32_t;
using word_t = std::uint64_t;
template <std::size_t K>
using edge_t = std::array<node_t, K>;
inline constexpr std::size_t WordBits = 64;

namespace detail
{

class Graph
{
public:
  explicit Graph(std::vector<node_t> vertices)
    : vertices_(std::move(vertices)), words_((vertices_.size() + WordBits - 1) / WordBits)
  {
    if(words_ && vertices_.size() > adjacency_.max_size() / words_) {
      throw std::length_error("PKMC adjacency matrix is too large");
    }
    adjacency_.assign(vertices_.size() * words_, 0);
  }

  std::size_t size() const noexcept
  {
    return vertices_.size();
  }

  node_t vertex(node_t local) const noexcept
  {
    return vertices_[local];
  }

  std::span<word_t> adjacency(node_t local) noexcept
  {
    return std::span(adjacency_).subspan(std::size_t(local) * words_, words_);
  }

  std::span<const word_t> adjacency(node_t local) const noexcept
  {
    return std::span(adjacency_).subspan(std::size_t(local) * words_, words_);
  }

  bool is_edge(node_t first, node_t second) const noexcept
  {
    return adjacency(first)[second / WordBits] & (word_t{1} << (second % WordBits));
  }

  void add_edge(node_t first, node_t second) noexcept
  {
    adjacency(first)[second / WordBits] |= word_t{1} << (second % WordBits);
    adjacency(second)[first / WordBits] |= word_t{1} << (first % WordBits);
  }

private:
  std::vector<node_t> vertices_;
  std::size_t words_;
  std::vector<word_t> adjacency_;
};

struct FaceHash
{
  template <std::size_t N>
  std::size_t operator()(const edge_t<N> &face) const noexcept
  {
    std::size_t hash = 0;
    for(node_t vertex : face) {
      hash ^= std::size_t(vertex) + 0x9e3779b9U + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

struct NoLinks
{
};

class Neighbors
{
public:
  void add(node_t vertex)
  {
    sparse_.push_back(vertex);
  }

  void pack(std::size_t vertex_count)
  {
    std::sort(sparse_.begin(), sparse_.end());
    const auto words = (vertex_count + WordBits - 1) / WordBits;
    if(sparse_.size() > words * sizeof(word_t) / sizeof(node_t)) {
      dense_.assign(words, 0);
      for(node_t vertex : sparse_) {
        dense_[vertex / WordBits] |= word_t{1} << (vertex % WordBits);
      }
      std::vector<node_t>().swap(sparse_);
    }
  }

  bool contains(node_t vertex) const noexcept
  {
    return dense_.empty() ? std::binary_search(sparse_.begin(), sparse_.end(), vertex)
                          : bool(dense_[vertex / WordBits] & (word_t{1} << (vertex % WordBits)));
  }

private:
  std::vector<node_t> sparse_;
  std::vector<word_t> dense_;
};

template <class Callback>
bool combinations(std::span<const node_t> vertices, std::span<node_t> face, Callback callback)
{
  const auto count = face.size();
  if(count > vertices.size()) {
    return true;
  }
  const auto visit = [&](auto &&self, std::size_t start, std::size_t depth) -> bool {
    if(depth == count) {
      return callback();
    }
    for(std::size_t i = start; i <= vertices.size() - (count - depth); ++i) {
      face[depth] = vertices[i];
      if(!self(self, i + 1, depth + 1)) {
        return false;
      }
    }
    return true;
  };
  return visit(visit, 0, 0);
}

}

template <std::size_t K>
class Hypergraph
{
  static_assert(K >= 2, "PKMC requires K >= 2");

public:
  using edge_type = edge_t<K>;
  using face_type = edge_t<K - 1>;

  Hypergraph(node_t vertex_count, std::span<const edge_type> edges)
    : edges_(edges.begin(), edges.end()), incidence_(vertex_count, 0),
      projection_(vertices(vertex_count))
  {
    for(auto &edge : edges_) {
      std::sort(edge.begin(), edge.end());
      if(edge.back() >= vertex_count) {
        throw std::out_of_range("PKMC vertex out of range");
      }
      if(std::adjacent_find(edge.begin(), edge.end()) != edge.end()) {
        throw std::invalid_argument("Repeated vertex in PKMC edge");
      }
    }
    std::sort(edges_.begin(), edges_.end());
    edges_.erase(std::unique(edges_.begin(), edges_.end()), edges_.end());
    for(const auto &edge : edges_) {
      for(std::size_t i = 0; i < K; ++i) {
        ++incidence_[edge[i]];
        for(std::size_t j = i + 1; j < K; ++j) {
          projection_.add_edge(edge[i], edge[j]);
        }
        if constexpr(K > 2) {
          face_type face{};
          for(std::size_t j = 0, index = 0; j < K; ++j) {
            if(j != i) {
              face[index++] = edge[j];
            }
          }
          links_[face].add(edge[i]);
        }
      }
    }
    if constexpr(K > 2) {
      for(auto &[key, neighbors] : links_) {
        neighbors.pack(vertex_count);
      }
    }
  }

  std::size_t size() const noexcept
  {
    return incidence_.size();
  }

  static constexpr std::size_t order() noexcept
  {
    return K;
  }

  const std::vector<edge_type> &edges() const noexcept
  {
    return edges_;
  }

  std::size_t incidence(node_t vertex) const noexcept
  {
    return incidence_[vertex];
  }

  const detail::Graph &projection() const noexcept
  {
    return projection_;
  }

  const detail::Neighbors *link(const face_type &sorted_face) const requires(K > 2)
  {
    const auto it = links_.find(sorted_face);
    return it == links_.end() ? nullptr : &it->second;
  }

  bool is_clique(std::span<const node_t> clique) const
  {
    std::vector<node_t> sorted(clique.begin(), clique.end());
    std::sort(sorted.begin(), sorted.end());
    if((!sorted.empty() && sorted.back() >= size()) ||
       std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
      return false;
    }
    edge_type edge{};
    return detail::combinations(sorted, edge, [&] {
      if constexpr(K == 2) {
        return projection_.is_edge(edge[0], edge[1]);
      } else {
        return std::binary_search(edges_.begin(), edges_.end(), edge);
      }
    });
  }

  // The caller supplies an existing clique and a distinct, in-range vertex.
  bool can_extend(std::span<const node_t> clique, node_t vertex) const
  {
    if constexpr(K == 2) {
      return std::ranges::all_of(
        clique, [&](node_t selected) { return projection_.is_edge(selected, vertex); });
    } else {
      face_type face{};
      return detail::combinations(clique, face, [&] {
        auto sorted_face = face;
        std::sort(sorted_face.begin(), sorted_face.end());
        const auto *neighbors = link(sorted_face);
        return neighbors && neighbors->contains(vertex);
      });
    }
  }

private:
  static std::vector<node_t> vertices(node_t count)
  {
    std::vector<node_t> result(count);
    std::iota(result.begin(), result.end(), node_t{0});
    return result;
  }

  std::vector<edge_type> edges_;
  std::vector<std::size_t> incidence_;
  detail::Graph projection_;
  [[no_unique_address]] std::conditional_t<K == 2,
    detail::NoLinks,
    std::unordered_map<face_type, detail::Neighbors, detail::FaceHash>>
    links_;
};

}
