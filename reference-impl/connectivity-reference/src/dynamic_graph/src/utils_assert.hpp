#pragma once

#include <cstdlib>
#include <iostream>
#include <sstream>

// If boolean `cond` is false, prints `msg` to stderr and terminates the
// program.
#define ASSERT_MSG_ALWAYS(cond, msg)                             \
  do {                                                           \
    if (!(cond)) {                                               \
      std::ostringstream str;                                    \
      str << msg;                                                \
      std::cerr << "ASSERTION FAILED: " << str.str() << '\n'     \
                << "at " << __FILE__ << ":" << __LINE__ << '\n'; \
      std::abort();                                              \
    }                                                            \
  } while (false)

#ifndef NDEBUG
#define ASSERT_MSG(cond, msg) ASSERT_MSG_ALWAYS(cond, msg)
#else
#define ASSERT_MSG(cond, msg) \
  do {                        \
  } while (false)
#endif  // ifndef NDEBUG
