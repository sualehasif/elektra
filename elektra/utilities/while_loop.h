#pragma once

#include <parlay/parallel.h>

namespace elektra {

// Launches several instances of function f, as if running a simple while loop
// (or do-while loop if is_do_while is true):
//     while (cond(state)) {
//       cilk_spawn f(state);
//       state = update_state(state);
//     }
//     cilk_sync;
template <typename Cond, typename F, typename UpdateState, typename State>
void ParallelWhile(Cond cond, F f, UpdateState update_state, State&& state, bool is_do_while = false) {
  if (is_do_while || cond(state)) {
    parlay::par_do(
      [&]() { f(state); },
      [&]() { ParallelWhile(cond, f, update_state, update_state(state)); }
    );
  }
}

}
