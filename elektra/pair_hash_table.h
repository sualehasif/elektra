#include <parlay/hash_table.h>

#include "hash_pair.h"

namespace elektra {

using std::pair;
using std::tuple;

template <typename T>
class PairHashTable {
 public:
  PairHashTable(int N) {
    table = new parlay::hashtable<parlay::hash_numeric<T>>(
        N, parlay::hash_numeric<T>{});
  }
  ~PairHashTable() { delete table; }

  // Inserts a new pair into the table.
  void insert(pair<T, T> p) { table->insert(hashIntPair(p)); }
  void insert(T p1, T p2) {
    table->insert(hashIntPair(std::make_pair(p1, p2)));
  }

  // Returns 1 if the table contains the pair.
  bool contains(pair<int, int> p) {
    auto hash = hashIntPair(p);
    return (table->find(hash) == hash);
  }

  void remove(pair<T, T> p) { table->deleteVal(hashIntPair(p)); }

 private:
  parlay::hashtable<parlay::hash_numeric<T>>* table;
};

// The struct for how to hash an edge
// K, V must be some integer type
// To be used with Parlay's hashtable
template <class K, class V>
struct hash_edge {
  using eType = tuple<K, K, V>;
  using kType = tuple<K, K, V>;
  eType empty() { return make_tuple(-1, -1, -1); }
  kType getKey(eType v) { return v; }
  size_t hash(kType v) {
    return static_cast<size_t>(
        hash64(std::get<0>(v) ^ std::get<1>(v) ^ std::get<2>(v)));
  }
  int cmp(kType v, kType b) {
    auto cmpSecond(kType v, kType b) =
        (std::get<1>(v) > std::get<1>(b))
            ? 1
            : ((std::get<1>(v) == std::get<1>(b)) ? 0 : -1);
    return (std::get<0>(v) > std::get<0>(b))
               ? 1
               : ((std::get<0>(v) == std::get<0>(b)) ? cmpSecond(v, b) : -1);
  }
  bool replaceQ(eType, eType) { return 0; }
  eType update(eType v, eType) { return v; }
  bool cas(eType* p, eType o, eType n) {
    // TODO: Make this use atomics properly. This is a quick
    // fix to get around the fact that the hashtable does
    // not use atomics. This will break for types that
    // do not inline perfectly inside a std::atomic (i.e.,
    // any type that the standard library chooses to lock)
    return std::atomic_compare_exchange_strong_explicit(
        reinterpret_cast<std::atomic<eType>*>(p), &o, n,
        std::memory_order_relaxed, std::memory_order_relaxed);
  }
};

}  // namespace elektra