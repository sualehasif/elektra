#pragma once

#include <limits.h>
#include <memory>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <fstream>   // std::ifstream
#include <iostream>  // std::cout
#include <tuple>

#include "utilities/cas.h"

namespace elektra {
//#define LONG 1

#ifndef NDEBUG
#define debug(_body) _body;
#else
#define debug(_body)
#endif

using uint = unsigned int;
using ulong = unsigned long;

// vertex size macros
using intV = int;
using uintV = unsigned int;

// size of edge-offsets.
// If the number of edges is more than sizeof(MAX_UINT),
// you should set the LONG flag on the command line.
#if defined(LONG)
using intT = long;
using uintT = unsigned long;
#define INT_T_MAX LONG_MAX
#define UINT_T_MAX ULONG_MAX
#else
using intT = int;
using uintT = unsigned int;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX
#endif

// edge size macros.
// If the number of vertices is more than sizeof(MAX_UINT)
// you should set the EDGELONG flag on the command line.
#if defined(EDGELONG)
typedef long intE;
typedef unsigned long uintE;
#define INT_E_MAX LONG_MAX
#define UINT_E_MAX ULONG_MAX
#else
using intE = int;
using uintE = unsigned int;
#define INT_E_MAX INT_MAX
#define UINT_E_MAX UINT_MAX
#endif

struct vertex_data {
  size_t offset;  // offset into the edges (static)
  uintE degree;   // possibly decreased by a (mutable) algorithm.
};

// Default granularity of a parallel for loop.
constexpr const size_t kDefaultGranularity = 2048;

#define LAST_BIT_SET(b) (b & (0x80))
#define EDGE_SIZE_PER_BYTE 7

struct empty {
  // equality operator
  auto operator==(const empty & /*unused*/) const -> bool { return true; }
};  // struct containing no data (used for empty base optimization)

// A templated sequence for ease of use.

using parlay::sequence;

}  // namespace elektra