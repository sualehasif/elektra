#pragma once

#include <parlay/parallel.h>
#include <parlay/utilities.h>

#include <tuple>

using namespace std;

template <class T>
struct maybe {
  T value;
  bool valid;

  maybe(T v, bool u) : value(v) { valid = u; }
  maybe(T v) : value(v) { valid = true; }
  maybe() { valid = false; }

  bool operator!() const { return !valid; }
  operator bool() const { return valid; };
  T &operator*() { return value; }
};

// #define granular_for(_i, _start, _end, _cond, _body)                 \
//   {                                                                  \
//     if (_cond) {                                                     \
//       {                                                              \
//         parlay::parallel_for(_start, _end, [&](size_t _i)) { _body } \
//       }                                                              \
//     } else {                                                         \
//       {                                                              \
//         for (size_t _i = _start; _i < _end; _i++) {                  \
//           _body                                                      \
//         }                                                            \
//       }                                                              \
//     }                                                                \
//   }

#define newA(__E, __n) (__E *)malloc((__n) * sizeof(__E))

#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX
#else
typedef int intT;
typedef unsigned int uintT;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX
#endif

template <class ET>
inline bool CAS(ET *ptr, ET oldv, ET newv) {
  if (sizeof(ET) == 1) {
    return __sync_bool_compare_and_swap((bool *)ptr, *((bool *)&oldv),
                                        *((bool *)&newv));
  } else if (sizeof(ET) == 4) {
    return __sync_bool_compare_and_swap((int *)ptr, *((int *)&oldv),
                                        *((int *)&newv));
  } else if (sizeof(ET) == 8) {
    return __sync_bool_compare_and_swap((long *)ptr, *((long *)&oldv),
                                        *((long *)&newv));
  }
#if defined(MCX16)
  else if (sizeof(ET) == 16) {
    return CAS128(ptr, oldv, newv);
  }
#endif
  else {
    std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
    abort();
  }
}

namespace concurrent_map {

// Thresholds for resizing down and up respectively
constexpr const double _load_lb = 0.05;
constexpr const double _load_ub = 0.95;
constexpr const double _growth_ratio = 2;
constexpr const size_t _minimum_ht_size = 4;

template <class K, class V, class HashFn>
class concurrentHT {
 public:
  using KV = tuple<K, V>;
  HashFn H;
  KV *table;
  K empty_key;
  K tombstone;
  bool alloc;

  size_t n_elms;
  size_t n_tombstones;
  size_t capacity;
  size_t mask;

  inline void clearA(KV *v, size_t n, K emp_key) {
    // granular_for(i, 0, n, (n > 2000), { get<0>(v[i]) = emp_key; });
    if (n > 2000) {
      parlay::parallel_for(0, n, [&](size_t i) { get<0>(v[i]) = emp_key; });
    } else {
      for (size_t i = 0; i < n; i++) {
        get<0>(v[i]) = emp_key;
      }
    }
  };

  inline size_t toRange(size_t h) const { return h & mask; }
  inline size_t firstIndex(K k) const { return toRange(H(k)); }
  inline size_t incrementIndex(size_t h) const { return toRange(h + 1); }

  inline KV *alloc_table(size_t _m) {
    // Must initialize for std::function()
    KV *tab = newA(KV, _m);
    clearA(tab, _m, empty_key);
    return tab;
  }

  inline void del() {
    if (alloc) {
      free(table);
    }
  }

  concurrentHT() : capacity(0), mask(0), table(0), n_elms(0), alloc(false) {}

  // Assumes size is a power of two.
  concurrentHT(KV *_table, size_t size, K _ek, K _ts)
      : table(_table),
        empty_key(_ek),
        tombstone(_ts),
        alloc(false),
        n_elms(0),
        n_tombstones(0),
        capacity(size),
        mask(size - 1) {
    if (table == nullptr) {
      capacity = 1 << parlay::log2_up(100 + (intT)(1.1 * (float)size));
      mask = capacity - 1;
      table = alloc_table(capacity);
      alloc = true;
    }
  }

  inline maybe<V> find(K k) const {
    size_t h = firstIndex(k);
    while (1) {
      KV t_kv = table[h];
      K t_k = get<0>(t_kv);
      if (t_k == k) {
        return maybe<V>(get<1>(t_kv));
      } else if (t_k == empty_key) {
        return maybe<V>();
      }
      h = incrementIndex(h);
    }
  }

  // Phase concurrent
  inline bool insert(K k, V v) {
    size_t h = firstIndex(k);
    while (1) {
      KV t_kv = table[h];
      K t_k = get<0>(t_kv);
      if ((t_k == empty_key && CAS(&(get<0>(table[h])), empty_key, k)) ||
          (t_k == tombstone && CAS(&(get<0>(table[h])), tombstone, k))) {
        get<1>(table[h]) = v;
        return true;
      } else if (t_k == k) {
        return false;
      }
      h = incrementIndex(h);
    }
  }

  // Phase concurrent
  inline bool deleteVal(K k) {
    size_t h = firstIndex(k);
    while (1) {
      KV t_kv = table[h];
      K t_k = get<0>(t_kv);
      if (t_k == empty_key) {
        return false;
      } else if (t_k == k) {
        // No atomics necessary?
        get<0>(table[h]) = tombstone;
        return true;
      }
      h = incrementIndex(h);
    }
  }
};

};  // namespace concurrent_map
