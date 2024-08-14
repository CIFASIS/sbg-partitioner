#include <fstream>
#include <iostream>


#pragma once

namespace sbg_partitioner {

class SBGLogger {
  public:
  static SBGLogger& instance()
  {
    static SBGLogger _instance;
    return _instance;
  }

  ~SBGLogger();

  std::ofstream log;

  private:
  SBGLogger();
};

static SBGLogger& sbg_log;

}