#pragma once

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/slice.h>

#include <tuple>
#include <unordered_set>
#include <vector>

#include "macros.h"

namespace elektra {

constexpr size_t kResizableTableCacheLineSz = 128;

inline size_t hashToRange(const size_t &h, const size_t &mask) {
  return h & mask;
}
inline size_t incrementIndex(const size_t &h, const size_t &mask) {
  return hashToRange(h + 1, mask);
}

template <class K, class V> struct iter_kv {
  using T = std::tuple<K, V>;
  K k;
  size_t h;
  size_t mask;
  // size_t num_probes;  // For debugging.
  T *table;
  K empty_key;
  iter_kv(K _k, size_t _h, size_t _mask, T *_table, K _empty_key)
      : k(_k), h(_h), mask(_mask),
        // num_probes(0),
        table(_table), empty_key(_empty_key) {}

  // Finds the location of the first key
  bool init() {
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        return false;
      } else if (std::get<0>(table[h]) == k) {
        return true;
      }
      h = incrementIndex(h, mask);
      // num_probes++;
    }
    return false;
  }

  // Probes until we've found another key
  bool has_next() {
    h = incrementIndex(h, mask);
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        return false;
      } else if (std::get<0>(table[h]) == k) {
        return true;
      }
      h = incrementIndex(h, mask);
      // num_probes++;
    }
    return false;
  }

  T next() { return table[h]; }
};

template <class K, class V, class KeyHash> class resizable_table {
public:
  using T = std::tuple<K, V>;

  size_t m;
  size_t mask;
  size_t ne;
  T empty;
  K empty_key;
  V empty_value;
  T tombstone;
  sequence<T> table;
  KeyHash key_hash;
  sequence<size_t> cts;

  inline size_t firstIndex(K &k) const { return hashToRange(key_hash(k), mask); }

  void init_counts() {
    size_t workers = parlay::num_workers();
    cts = sequence<size_t>::uninitialized(kResizableTableCacheLineSz * workers);
    for (size_t i = 0; i < workers; i++) {
      cts[i * kResizableTableCacheLineSz] = 0;
    }
  }

  void update_nelms() {
    size_t workers = parlay::num_workers();
    for (size_t i = 0; i < workers; i++) {
      ne += cts[i * kResizableTableCacheLineSz];
      cts[i * kResizableTableCacheLineSz] = 0;
    }
  }

  resizable_table() : m(0), ne(0) {
    mask = 0;
    init_counts();
  }

  resizable_table(size_t _m, T _empty, T _tombstone, KeyHash _key_hash)
      : m((size_t)1 << parlay::log2_up((size_t)(1.1 * _m))), mask(m - 1), ne(0),
        empty(_empty), empty_key(std::get<0>(empty)),
        empty_value(std::get<1>(empty)), tombstone(_tombstone),
        key_hash(_key_hash) {
    table = sequence<T>::uninitialized(m);
    clear();
    init_counts();
  }

  void analyze() {
    std::vector<size_t> probe_lengths;
    size_t chain_length = 0;
    for (size_t i = 0; i < m; i++) {
      if (std::get<0>(table[i]) == empty_key ||
          std::get<0>(table[i]) == tombstone) {
        if (chain_length > 0) {
          probe_lengths.emplace_back(chain_length);
        }
        chain_length = 0;
      } else {
        chain_length++;
      }
    }

    std::sort(probe_lengths.begin(), probe_lengths.end());
    std::cout << "Analyzed table, m = " << m << " ne = " << ne
              << " num_probes = " << probe_lengths.size() << std::endl;
    for (long k = probe_lengths.size() - 1; k > 0; k--) {
      std::cout << probe_lengths[k] << " ";
      if (probe_lengths.size() - k > 250)
        break;
    }
    std::cout << std::endl << "End of table analysis" << std::endl;
  }

  void maybe_resize(size_t n_inc) {
    update_nelms();
    size_t nt = ne + n_inc;
    if (nt > (0.25 * m)) {
      size_t old_m = m;
      auto old_t = std::move(table);
      m = ((size_t)1 << parlay::log2_up((size_t)(10 * nt)));
      if (m == old_m) {
        return;
      }
      mask = m - 1;
      ne = 0;
      table = sequence<T>::uninitialized(m);
      clear();
      parlay::parallel_for(0, old_m, [&](size_t i) {
        if (std::get<0>(old_t[i]) != empty_key &&
            std::get<0>(old_t[i]) != std::get<0>(tombstone)) {
          insert(old_t[i]);
        }
      });
      update_nelms();
    }
  }

  iter_kv<K, V> get_iter(K k) {
    size_t h = firstIndex(k);
    return iter_kv<K, V>(k, h, mask, table.begin(), empty_key);
  }

  // TODO URGENT what about tombstones?
  bool insert(std::tuple<K, V> kv) {
    K &k = std::get<0>(kv);
    V &v = std::get<1>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == empty_key &&
          elektra::atomic_compare_and_swap(&std::get<0>(table[h]), empty_key,
                                           std::get<0>(kv))) {
        std::get<1>(table[h]) = v;
        size_t wn = parlay::worker_id();
        cts[wn * kResizableTableCacheLineSz]++;
        return 1;
      }
      if (std::get<0>(table[h]) == k) {
        std::get<1>(table[h]) = v;
        return 1;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  bool insert_seq(std::tuple<K, V> kv) {
    K &k = std::get<0>(kv);
    V &v = std::get<1>(kv);
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == empty_key) {
        table[h] = kv;
        return 1;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  size_t num_appearances(K k) {
    size_t h = firstIndex(k);
    size_t ct = 0;
    while (1) {
      if (std::get<0>(table[h]) == k) {
        ct++;
      } else if (std::get<0>(table[h]) == empty_key) {
        return ct;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  V find(K k) {
    size_t h = firstIndex(k);

    while (1) {
      if (std::get<0>(table[h]) == k) {
        return std::get<1>(table[h]);
      } else if (std::get<0>(table[h]) == empty_key) {
        return empty_value;
      }
      h = incrementIndex(h, mask);
    }
    return empty_value;
  }

  bool contains(K k) const {
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == k) {
        return 1;
      } else if (std::get<0>(table[h]) == empty_key) {
        return 0;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  bool contains(K k, V v) const {
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == k) {
        return std::get<1>(table[h]) == v;
      } else if (std::get<0>(table[h]) == empty_key) {
        return 0;
      }
      h = incrementIndex(h, mask);
    }
    return 0;
  }

  void deleteVal(K k) {
    size_t h = firstIndex(k);
    while (1) {
      if (std::get<0>(table[h]) == k) {
        table[h] = tombstone;
        return;
      } else if (std::get<0>(table[h]) == empty_key) {
        return;
      }
      h = incrementIndex(h, mask);
    }
  }

  template <class F> void map(F &f) {
    parallel_for(0, m, [&](size_t i) {
      if (std::get<0>(table[i]) != empty_key &&
          std::get<0>(table[i]) != std::get<0>(tombstone)) {
        f(table[i]);
      }
    });
  }

  sequence<T> entries() const {
    auto pred = [&](const T &t) {
      return (std::get<0>(t) != empty_key &&
              std::get<0>(t) != std::get<0>(tombstone));
    };
    auto table_seq = parlay::make_slice(table);
    return parlay::filter(table_seq, pred);
  }

  sequence<K> keys() const {
    auto pred = [&](const K &k) {
      return k != empty_key && k != std::get<0>(tombstone);
    };
    auto table_seq = parlay::make_slice(table);
    auto keys_seq = parlay::delayed_seq<K>(table_seq.size(), [&](const size_t i) {
        return std::get<0>(table_seq[i]);
    });
    return parlay::filter(keys_seq, pred);
  }

  void clear() {
    parlay::parallel_for(0, m, [&](size_t i) { table[i] = empty; });
  }

  size_t size() {
    update_nelms();
    return ne;
  }
};

template <class K, class V, class KeyHash>
inline resizable_table<K, V, KeyHash>
make_resizable_table(size_t m, std::tuple<K, V> empty,
                     std::tuple<K, V> tombstone, KeyHash key_hash) {
  return resizable_table<K, V, KeyHash>(m, empty, tombstone, key_hash);
}
} // namespace elektra
