// Concurrent array allocator for small arrays.
// Elements of the array are not initialized. If T is an object, the allocator
// does not call constructors or destructors for T.
//
// Do not statically initialize this. This object depends on list_allocator
// being statically initialized before it can be initialized.
//
// Implementation: to handle an allocation request of length k, round k up to
// the nearest power of 2 and give an array of that size. Thus if the max array
// size is n, there are log(n) sizes to allocate. We handle each of these sizes
// with a concurrent fixed-size allocator.
#pragma once

#include <parlay/alloc.h>
#include <parlay/utilities.h>

#include <cassert>

namespace concurrent_array_allocator {

constexpr int kMaxArrayLength{32};

template <typename T>
class Allocator {
 public:
  Allocator();
  // Note that this destructor will call `finish()` on several
  // `list_allocator<T>`s, so be sure that allocators of the same type aren't
  // used elsewhere.
  // TODO(tomtseng): We shouldn't have to think about other objects when calling
  // a destructor.
  ~Allocator();
  Allocator(const Allocator&) = delete;
  Allocator(Allocator&&) = delete;
  Allocator& operator=(const Allocator&) = delete;
  Allocator& operator=(Allocator&&) = delete;

  T* Allocate(int length);
  void Free(T* arr, int length);

 private:
  using list_allocator0 = parlay::type_allocator<T[1]>;
  using list_allocator1 = parlay::type_allocator<T[2]>;
  using list_allocator2 = parlay::type_allocator<T[4]>;
  using list_allocator3 = parlay::type_allocator<T[8]>;
  using list_allocator4 = parlay::type_allocator<T[16]>;
  using list_allocator5 = parlay::type_allocator<T[32]>;
};

///////////////////////////////////////////////////////////////////////////////
//                           Implementation below.                           //
///////////////////////////////////////////////////////////////////////////////

template <typename T>
Allocator<T>::Allocator() {}

template <typename T>
Allocator<T>::~Allocator() {}

template <typename T>
T* Allocator<T>::Allocate(int length) {
  switch (parlay::log2_up(length)) {
    case 0:
      return *list_allocator0::alloc();
    case 1:
      return *list_allocator1::alloc();
    case 2:
      return *list_allocator2::alloc();
    case 3:
      return *list_allocator3::alloc();
    case 4:
      return *list_allocator4::alloc();
    case 5:
      return *list_allocator5::alloc();
    default:
      return nullptr;
  }
}

template <typename T>
void Allocator<T>::Free(T* arr, int length) {
  // To justify this `reinterpret_cast`, we use that a pointer to a static array
  // and a static array itself both point to the same address: the base of the
  // array.
  switch (parlay::log2_up(length)) {
    case 0:
      list_allocator0::free(reinterpret_cast<T(*)[1]>(arr));
      break;
    case 1:
      list_allocator1::free(reinterpret_cast<T(*)[2]>(arr));
      break;
    case 2:
      list_allocator2::free(reinterpret_cast<T(*)[4]>(arr));
      break;
    case 3:
      list_allocator3::free(reinterpret_cast<T(*)[8]>(arr));
      break;
    case 4:
      list_allocator4::free(reinterpret_cast<T(*)[16]>(arr));
      break;
    case 5:
      list_allocator5::free(reinterpret_cast<T(*)[32]>(arr));
      break;
    default:
      assert(false);
  }
}

}  // namespace concurrent_array_allocator