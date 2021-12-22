#pragma once

#include <parlay/sequence.h>

#include "macros.h"

namespace elektra {

/**************************** Framework options ****************************/
// enum SamplingOption { sample_kout, sample_bfs, sample_ldd, no_sampling };

/* Union-Find options */
enum FindOption {
  find_compress,
  find_naive,
  find_split,
  find_halve,
  find_atomic_split,
  find_atomic_halve
};
enum UniteOption {
  unite,
  unite_early,
  unite_nd,
  unite_rem_lock,
  unite_rem_cas
};

/* RemCAS-specific options */
enum SpliceOption {
  split_atomic_one,
  halve_atomic_one,
  splice_simple,
  splice_atomic
};

enum UpdateType { insertion_type, query_type };

using parent = uintE;

void report_pathlen(uintE pathlen) {
#ifdef REPORT_PATH_LENGTHS
  max_pathlen.update_value(pathlen);
  total_pathlen.update_value(pathlen);
#endif
}

namespace jayanti_rank {
static constexpr uintE RANK_MASK = (uintE) INT_E_MAX;
static constexpr uintE TOP_BIT_SHIFT = sizeof(uintE) * 8 - 1;
static constexpr uintE TOP_BIT = ((uintE) 1) << TOP_BIT_SHIFT;

struct vdata {
  volatile uintE rank;  // top bit is used to indicate root or not
  volatile uintE par;   // parent id

  vdata() {}

  vdata(uintE _parent, uintE _rank, bool is_root) {
    rank = combine_root_rank(is_root, _rank);
    par = _parent;
  }

  __attribute__((always_inline)) inline uintE combine_root_rank(
      bool is_root, uintE _rank) const {
    return (((uintE) is_root) << TOP_BIT_SHIFT) + _rank;
  }

  __attribute__((always_inline)) inline bool is_root() const {
    return rank & TOP_BIT;
  }

  __attribute__((always_inline)) inline uintE get_rank() const {
    return rank & RANK_MASK;
  }

  __attribute__((always_inline)) inline uintE get_parent() const { return par; }

  void print(uintE vtx_id) const {
    std::cout << "# vtx: " << vtx_id << " parent = " << get_parent()
              << " rank = " << get_rank() << " is_root = " << is_root()
              << std::endl;
  }
};

template<class S>
void link(uintE u, uintE v, S &vdatas, parlay::random r) {
  auto ud = vdatas[u];
  auto vd = vdatas[v];
  // spend two reads to abort early with no CASs if either of the
  // vertices has already been updated.
  if (ud.is_root() == false || vd.is_root() == false) {
    return;
  }
  // o.w. continue
  if (ud.get_rank() < vd.get_rank()) {  // u.r < v.r
    auto expected_u =
        vdata(ud.get_parent(), ud.get_rank(), /* is_root= */ true);
    auto new_u = vdata(v, ud.get_rank(), /* is_root= */ false);
    atomic_compare_and_swap(&(vdatas[u]), expected_u, new_u);
  } else if (ud.get_rank() > vd.get_rank()) {  // v.r < u.r
    auto expected_v =
        vdata(vd.get_parent(), vd.get_rank(), /* is_root= */ true);
    auto new_v = vdata(u, vd.get_rank(), /* is_root= */ false);
    atomic_compare_and_swap(&(vdatas[v]), expected_v, new_v);
  } else {  // u.r == v.r
    auto random_bit = r.rand() & 1;
    if (u < v) {
      auto expected_u =
          vdata(ud.get_parent(), ud.get_rank(), /* is_root= */ true);
      auto new_u = vdata(v, ud.get_rank() + 1, /* is_root= */ random_bit);
      atomic_compare_and_swap(&(vdatas[u]), expected_u, new_u);
    } else {
      auto expected_v =
          vdata(vd.get_parent(), vd.get_rank(), /* is_root= */ true);
      auto new_v = vdata(u, vd.get_rank() + 1, /* is_root= */ random_bit);
      atomic_compare_and_swap(&(vdatas[v]), expected_v, new_v);
    }
  }
}

inline uintE find(uintE x, sequence<vdata> &vdatas) {
  uintE u = x;
  uintE pathlen = 1;
  while (!vdatas[u].is_root()) {  // * on u.is_root()
    u = vdatas[u].get_parent();
    pathlen++;
  }
  report_pathlen(pathlen);
  return u;  // u is a root
}

inline uintE find_twotry_splitting(uintE x, sequence<vdata> &vdatas) {
  uintE u = x;
  uintE pathlen = 1;
  while (!vdatas[u].is_root()) {  // * on u.is_root()
    auto ud = vdatas[u];
    uintE v = ud.get_parent();
    auto vd = vdatas[v];
    pathlen++;
    if (vd.is_root()) {
      report_pathlen(pathlen);
      return v;
    }

    // CAS 1
    uintE w = vd.get_parent();
    auto expected_u = vdata(v, ud.get_rank(), false);
    auto new_u = vdata(w, ud.get_rank(), false);
    atomic_compare_and_swap<vdata>(&(vdatas[u]), expected_u, new_u);

    // read and check
    ud = vdatas[u];
    v = ud.get_parent();
    vd = vdatas[v];
    w = vd.get_parent();
    if (vd.is_root()) {
      report_pathlen(pathlen);
      return v;
    }

    pathlen++;

    // CAS 2
    expected_u = vdata(v, ud.get_rank(), false);
    new_u = vdata(w, ud.get_rank(), false);
    atomic_compare_and_swap<vdata>(&(vdatas[u]), expected_u, new_u);

    u = v;
  }
  report_pathlen(pathlen);
  return u;  // u is a root
}

template<class S, class Find>
void unite(uintE x, uintE y, S &vdatas, parlay::random r, Find &find) {
  uintE u = find(x, vdatas);
  uintE v = find(y, vdatas);
  while (u != v) {
    link(u, v, vdatas, r);
    u = find(u, vdatas);
    v = find(v, vdatas);
    r = r.next();
  }
}
};  // namespace jayanti_rank

namespace find_variants {
inline uintE find_naive(uintE i, sequence<parent> &parents) {
  while (i != parents[i]) {
    i = parents[i];
  }
  return i;
}

inline uintE find_compress(uintE i, sequence<parent> &parents) {
  parent j = i;
  if (parents[j] == j) return j;
  do {
    j = parents[j];
  } while (parents[j] != j);
  parent tmp;
  while ((tmp = parents[i]) > j) {
    parents[i] = j;
    i = tmp;
  }
  return j;
}

inline uintE find_atomic_split(uintE i, sequence<parent> &parents) {
  while (1) {
    parent v = parents[i];
    parent w = parents[v];
    if (v == w) {
      return v;
    } else {
      atomic_compare_and_swap(&parents[i], v, w);
      // i = its parents_
      i = v;
    }
  }
}

inline uintE find_atomic_halve(uintE i, sequence<parent> &parents) {
  while (1) {
    parent v = parents[i];
    parent w = parents[v];
    if (v == w) {
      return v;
    } else {
      atomic_compare_and_swap(&parents[i], (parent) v, (parent) w);
      // i = its grandparent
      i = parents[i];
    }
  }
}
}  // namespace find_variants

namespace splice_variants {

/* Used in Rem-CAS variants for splice */
inline uintE split_atomic_one(uintE i, uintE x, sequence<parent> &parents) {
  parent v = parents[i];
  parent w = parents[v];
  if (v == w)
    return v;
  else {
    atomic_compare_and_swap(&parents[i], v, w);
    i = v;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
inline uintE halve_atomic_one(uintE i, uintE x, sequence<parent> &parents) {
  parent v = parents[i];
  parent w = parents[v];
  if (v == w)
    return v;
  else {
    atomic_compare_and_swap(&parents[i], v, w);
    i = w;
    return i;
  }
}

/* Used in Rem-CAS variants for splice */
inline uintE splice_atomic(uintE u, uintE v, sequence<parent> &parents) {
  parent z = parents[u];
  atomic_compare_and_swap(&parents[u], z, parents[v]);
  return z;
}
}  // namespace splice_variants

namespace unite_variants {

template<class Find>
struct Unite {
  Find &find;
  Unite(Find &find) : find(find) {}

  using edge = std::pair<uintE, uintE>;

  inline uintE operator()(uintE u_orig, uintE v_orig,
                          sequence<parent> &parents) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, parents);
      v = find(v, parents);
      if (u == v)
        break;
      else if (u > v && parents[u] == u &&
          atomic_compare_and_swap(&parents[u], u, v)) {
        return u;
      } else if (v > u && parents[v] == v &&
          atomic_compare_and_swap(&parents[v], v, u)) {
        return v;
      }
    }
    return UINT_E_MAX;
  }

  inline void operator()(uintE u_orig, uintE v_orig, sequence<parent> &parents,
                         sequence<edge> &edges) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, parents);
      v = find(v, parents);
      if (u == v)
        break;
      else if (u > v && parents[u] == u &&
          atomic_compare_and_swap(&parents[u], u, v)) {
        edges[u] = std::make_pair(u_orig, v_orig);
        break;
      } else if (v > u && parents[v] == v &&
          atomic_compare_and_swap(&parents[v], v, u)) {
        edges[v] = std::make_pair(u_orig, v_orig);
        break;
      }
    }
  }
};

template<class Splice, class Compress, FindOption find_option>
struct UniteRemLock {
  uintE n;
  sequence<bool> locks;
  Compress &compress;
  Splice &splice;
  UniteRemLock(Compress &compress, Splice &splice, uintE n)
      : n(n), compress(compress), splice(splice) {
    locks = sequence<bool>(n, false);
  }

  inline void fence() { std::atomic_thread_fence(std::memory_order_seq_cst); }

  bool acquire_lock(uintE u) {
    if (!locks[u] && atomic_compare_and_swap(&locks[u], false, true)) {
      return true;
    }
    return false;
  }

  void release_lock(uintE u) { locks[u] = false; }

  inline uintE operator()(uintE u_orig, uintE v_orig,
                          sequence<parent> &parents) {
    parent rx = u_orig;
    parent ry = v_orig;
    while (parents[rx] != parents[ry]) {
      /* link from high -> low */
      if (parents[rx] < parents[ry]) std::swap(rx, ry);
      if (rx == parents[rx]) {
        if (acquire_lock(rx)) {
          parent py = parents[ry];
          if (rx == parents[rx] && rx > py) {
            parents[rx] = py;
          }
          release_lock(rx);
        }
        /* link high -> low */
      } else {
        rx = splice(rx, ry, parents);
      }
    }
    if constexpr (find_option != find_naive) { /* aka find_none */
      compress(u_orig, parents);
      compress(v_orig, parents);
    }
    return UINT_E_MAX;  // TODO (ret-value)
  }
};

template<class Splice, class Compress, FindOption find_option>
struct UniteRemCAS {
  Compress &compress;
  Splice &splice;
  UniteRemCAS(Compress &compress, Splice &splice)
      : compress(compress), splice(splice) {}

  inline uintE operator()(uintE x, uintE y, sequence<parent> &parents) {
    uintE rx = x;
    uintE ry = y;
    while (parents[rx] != parents[ry]) {
      /* link high -> low */
      parent p_ry = parents[ry];
      parent p_rx = parents[rx];
      if (p_rx < p_ry) {
        std::swap(rx, ry);
        std::swap(p_rx, p_ry);
      }
      if (rx == parents[rx] &&
          atomic_compare_and_swap(&parents[rx], rx, p_ry)) {
        if constexpr (find_option != find_naive) { /* aka find_none */
          compress(x, parents);
          compress(y, parents);
        }
        return rx;
      } else {
        // failure: locally compress by splicing and try again
        rx = splice(rx, ry, parents);
      }
    }
    return UINT_E_MAX;
  }
};

template<class Find, FindOption find_option>
struct UniteEarly {
  Find &find;
  UniteEarly(Find &find) : find(find) {}
  inline uintE operator()(uintE u, uintE v, sequence<parent> &parents) {
    [[maybe_unused]] uintE u_orig = u, v_orig = v;
    uintE ret = UINT_E_MAX;
    while (u != v) {
      /* link high -> low */
      if (v > u) std::swap(u, v);
      if (parents[u] == u && atomic_compare_and_swap(&parents[u], u, v)) {
        ret = u;
        break;
      }
      parent z = parents[u];
      parent w = parents[z];
      atomic_compare_and_swap(&parents[u], z, w);
      u = w;
    }
    if constexpr (find_option != find_naive) {
      u = find(u_orig, parents); /* force */
      v = find(v_orig, parents); /* force */
    }
    return ret;
  }
};

template<class Find>
struct UniteND {
  Find find;
  sequence<uintE> hooks;
  UniteND(size_t n, Find &find) : find(find) {
    hooks = sequence<uintE>(n, UINT_E_MAX);
  }

  inline uintE operator()(uintE u_orig, uintE v_orig,
                          sequence<parent> &parents) {
    parent u = u_orig;
    parent v = v_orig;
    while (1) {
      u = find(u, parents);
      v = find(v, parents);
      if (u == v) break;
      /* link high -> low */
      if (u < v) std::swap(u, v);
      if (hooks[u] == UINT_E_MAX &&
          atomic_compare_and_swap(&hooks[u], UINT_E_MAX, v)) {
        parents[u] = v;
        return u;
      }
    }
    return UINT_E_MAX;
  }
};

}  // namespace unite_variants

struct UnionFindCompress {
  unsigned long size;
  static constexpr auto find{find_variants::find_compress};
  unite_variants::Unite<decltype(find)> unite{unite_variants::Unite<decltype(find)>{find}};
  parlay::sequence<parent> parents = parlay::sequence<parent>::from_function(size, [](uintE i) { return i; });
};

}  // namespace elektra