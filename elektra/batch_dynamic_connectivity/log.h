#define EXPAND_LABEL(x) (#x)
#define EXPAND_VALUE(x) x

#define LOG_VAL(msg, label, value)                            \
  {                                                           \
    std::cout << std::string("[") + std::string(" ");         \
    std::cout << __FILE__;                                    \
    std::cout << std::string(":");                            \
    std::cout << __LINE__;                                    \
    std::cout << std::string("] ") + msg + std::string("\n"); \
    std::cout << EXPAND_LABEL(label);                         \
    std::cout.flush();                                        \
  }