#define PARENS ()

#define EXPAND(...) EXPAND4(EXPAND4(EXPAND4(EXPAND4(__VA_ARGS__))))
#define EXPAND4(...) EXPAND3(EXPAND3(EXPAND3(EXPAND3(__VA_ARGS__))))
#define EXPAND3(...) EXPAND2(EXPAND2(EXPAND2(EXPAND2(__VA_ARGS__))))
#define EXPAND2(...) EXPAND1(EXPAND1(EXPAND1(EXPAND1(__VA_ARGS__))))
#define EXPAND1(...) __VA_ARGS__

#define FOR_EACH(macro, between, ...) \
  __VA_OPT__(EXPAND(FOR_EACH_HELPER(macro, __VA_ARGS__)))
#define FOR_EACH_HELPER(macro, a1, ...) \
  macro(a1) __VA_OPT__(, FOR_EACH_AGAIN PARENS(macro, __VA_ARGS__))
#define FOR_EACH_AGAIN() FOR_EACH_HELPER

#define FOR_EACH2(macro1, macro2, ...) \
  __VA_OPT__(EXPAND(FOR_EACH2_HELPER(macro1, macro2, __VA_ARGS__)))
#define FOR_EACH2_HELPER(macro1, macro2, a1, a2, ...) \
  macro1(a1),                                         \
      macro2(a2)                                      \
          __VA_OPT__(, FOR_EACH2_AGAIN PARENS(macro1, macro2, __VA_ARGS__))
#define FOR_EACH2_AGAIN() FOR_EACH2_HELPER

#define EXPAND_LABEL(x) (" " + #x + "=")
#define EXPAND_VALUE(x) x
#define LOG(msg, ...)                                               \
  {                                                                           \
    auto s = ("[" + " " + __FILE__ + ":" \
                        + __LINE__ + "] "+ msg,                                \
                          __VA_OPT__(FOR_EACH2(EXPAND_LABEL, EXPAND_VALUE,    \
                                               __VA_ARGS__), ) "\n");         \
    std::cerr << s;                                                           \
    std::cerr.flush();                                                        \
  }