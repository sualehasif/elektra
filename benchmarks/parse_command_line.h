
#pragma once

#include <stdlib.h>

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace elektra {

struct commandLine {
  int argc;
  char** argv;
  std::string comLine;

  commandLine(int _c, char** _v, std::string _cl)
      : argc(_c), argv(_v), comLine(std::move(_cl)) {}

  commandLine(int _c, char** _v)
      : argc(_c), argv(_v), comLine("bad arguments") {}

  void badArgument() const {
    std::cout << "usage: " << argv[0] << " " << comLine << "\n";
    abort();
  }

  // get an argument
  // i is indexed from the last argument = 0, second to last indexed 1, ..
  char* getArgument(int i) const {
    if (argc < 2 + i) badArgument();
    return argv[argc - 1 - i];
  }

  // looks for two filenames
  std::pair<char*, char*> IOFileNames() const {
    if (argc < 3) badArgument();
    return std::pair<char*, char*>(argv[argc - 2], argv[argc - 1]);
  }

  std::pair<int, char*> sizeAndFileName() const {
    if (argc < 3) badArgument();
    return std::pair<int, char*>(atoi(argv[argc - 2]), (char*)argv[argc - 1]);
  }

  bool getOption(const std::string& option) const {
    for (int i = 1; i < argc; i++)
      if ((std::string)argv[i] == option) return true;
    return false;
  }

  char* getOptionValue(const std::string& option) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) return argv[i + 1];
    return NULL;
  }

  std::string getOptionValue(const std::string& option,
                             const std::string& defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) return (std::string)argv[i + 1];
    return defaultValue;
  }

  int getOptionIntValue(const std::string& option, int defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) {
        int r = atoi(argv[i + 1]);
        return r;
      }
    return defaultValue;
  }

  size_t getOptionLongValue(const std::string& option,
                            size_t defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) {
        long r = atol(argv[i + 1]);
        return r;
      }
    return defaultValue;
  }

  double getOptionDoubleValue(const std::string& option,
                              double defaultValue) const {
    for (int i = 1; i < argc - 1; i++)
      if ((std::string)argv[i] == option) {
        double val;
        if (sscanf(argv[i + 1], "%lf", &val) == EOF) {
          badArgument();
        }
        return val;
      }
    return defaultValue;
  }
};

}  // namespace elektra
