#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#ifdef LATTICE_CORE_TEST
#define LOGDBG(X)                                                              \
  {                                                                            \
    X                                                                          \
  }
#else
#define LOGDBG(X)
#endif

// #ifdef LATTICE_CORE_TEST
#define FLOGOUT(X)                                                             \
  {                                                                            \
    DebugTools::Log::get_instance().flog X << std::endl;                       \
  }
// #else
// #define FLOGOUT(X)
// #endif

namespace DebugTools {
class Log
{
private:
  std::string LOG_PATH;

  std::string flog_filename;

  const bool log_on = true;

  decltype(std::chrono::steady_clock::now()) time_p[2];

  int time_id = 0;

  Log();

  Log(const Log& log) = delete;
  Log(Log&& log) = delete;

protected:
public:
  ~Log() = default;

  std::ofstream flog;
  int file_count = 0;

  void create_new_log(std::string filename = "");

  void tick();

  double get_millisec_span();

  void print_time_span_tostd();

  static Log& get_instance()
  {
    static Log _log;
    return _log;
  }
};

} // namespace DebugTools