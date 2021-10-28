#include <parlay/hash_table.h>

#include "hash_pair.h"

namespace elektra {

using std::pair;

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
  parlay::hashtable<parlay::hash_numeric<T>> *table;
};
}  // namespace elektra