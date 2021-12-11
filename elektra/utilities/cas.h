#pragma once

#include <type_traits>
#include <utility>

template <class...> constexpr std::false_type always_false{};

namespace elektra {

// ========================= atomic ops  ==========================
template <typename ET>
inline bool atomic_compare_and_swap(ET* a, ET oldval, ET newval) {
  if constexpr (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint8_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(a), r_oval,
                                        r_nval);
  } else if constexpr (sizeof(ET) == 16) {
    __int128 r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap_16(reinterpret_cast<__int128*>(a),
                                           r_oval, r_nval);
  } else {
    static_assert(always_false<ET>);
  }
}

template <typename ET>
inline bool atomic_compare_and_swap(volatile ET* a, ET oldval, ET newval) {
  if constexpr (sizeof(ET) == 1) {
    uint8_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint8_t*>(a),
                                        r_oval, r_nval);
  } else if constexpr (sizeof(ET) == 4) {
    uint32_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint32_t*>(a),
                                        r_oval, r_nval);
  } else if constexpr (sizeof(ET) == 8) {
    uint64_t r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap(reinterpret_cast<volatile uint64_t*>(a),
                                        r_oval, r_nval);
  } else if constexpr (sizeof(ET) == 16) {
    __int128 r_oval, r_nval;
    std::memcpy(&r_oval, &oldval, sizeof(ET));
    std::memcpy(&r_nval, &newval, sizeof(ET));
    return __sync_bool_compare_and_swap_16(
        reinterpret_cast<volatile __int128*>(a), r_oval, r_nval);
  } else {
    static_assert(always_false<ET>);
  }
}

template <typename E, typename EV>
inline E fetch_and_add(E* a, EV b) {
  volatile E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
  return oldV;
}

template <typename E, typename EV>
inline void write_add(E* a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = *a;
    newV = oldV + b;
  } while (!atomic_compare_and_swap(a, oldV, newV));
}

template <typename E, typename EV>
inline void write_add(std::atomic<E>* a, EV b) {
  // volatile E newV, oldV;
  E newV, oldV;
  do {
    oldV = a->load();
    newV = oldV + b;
  } while (!std::atomic_compare_exchange_strong(a, &oldV, newV));
}

template <typename ET, typename F>
inline bool write_min(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_min(volatile ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = *a;
  while (less(b, c) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_min(std::atomic<ET>* a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = a->load();
  while (less(b, c) && !(r = std::atomic_compare_exchange_strong(a, &c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(volatile ET* a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = *a;
  while (less(c, b) && !(r = atomic_compare_and_swap(a, c, b)));
  return r;
}

template <typename ET, typename F>
inline bool write_max(std::atomic<ET>* a, ET b, F less) {
  ET c;
  bool r = 0;
  do c = a->load();
  while (less(c, b) && !(r = std::atomic_compare_exchange_strong(a, &c, b)));
  return r;
}

}  // namespace elektra

template <class ET>
inline bool CAS(ET* ptr, ET oldv, ET newv) {
  if constexpr (sizeof(ET) == 1) {
    return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv),
                                        *((bool*)&newv));
  } else if constexpr (sizeof(ET) == 4) {
    return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv),
                                        *((int*)&newv));
  } else if constexpr (sizeof(ET) == 8) {
    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv),
                                        *((long*)&newv));
  }
#if defined(MCX16)
  else if constexpr (sizeof(ET) == 16) {
    return CAS128(ptr, oldv, newv);
  }
#endif
  else {
    static_assert(always_false<ET>);
  }
}
