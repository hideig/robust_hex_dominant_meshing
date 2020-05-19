#include "log.h"

namespace DebugTools {
Log::Log()
{
#ifdef _WIN32
  LOG_PATH = "./";
#endif
#ifdef __linux
  LOG_PATH = "/dev/shm/";
#endif
  flog_filename = LOG_PATH + "/flog1.log";
  flog = std::ofstream(flog_filename);
}

void
Log::create_new_log(std::string filename)
{
  if (filename == "") {
    filename = LOG_PATH + "log" + std::to_string(file_count);
    file_count++;
  } else {
    filename = LOG_PATH + filename;
  }
  flog.close();
  flog.open(filename);
}

void
Log::tick()
{
  time_p[time_id] = std::chrono::steady_clock::now();
  time_id ^= 1;
}

double
Log::get_millisec_span()
{
  return std::chrono::duration<double, std::milli>(time_p[time_id ^ 1] -
                                                   time_p[time_id])
    .count();
}

void
Log::print_time_span_tostd()
{
  double millisecs = get_millisec_span();
  if (millisecs < 1000)
    std::cout << "time span: " << millisecs << " ms" << std::endl;
  else {
    millisecs /= 1000;
    std::cout << "time span: " << millisecs << " s" << std::endl;
  }
}
} // namespace DebugTools