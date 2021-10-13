#include <stdlib.h>
#include <sys/time.h>

#include <iomanip>
#include <iostream>

namespace elektra {
struct timer {
  double totalTime;
  double lastTime;
  double totalWeight;
  bool on;
  struct timezone tzp;
  timer() {
    struct timezone tz = {0, 0};
    totalTime = 0.0;
    totalWeight = 0.0;
    on = 0;
    tzp = tz;
  }
  double getTime() {
    timeval now;
    gettimeofday(&now, &tzp);
    return ((double)now.tv_sec) + ((double)now.tv_usec) / 1000000.;
  }
  void start() {
    on = 1;
    lastTime = getTime();
  }
  double stop() {
    on = 0;
    double d = (getTime() - lastTime);
    totalTime += d;
    return d;
  }
  double stop(double weight) {
    on = 0;
    totalWeight += weight;
    double d = (getTime() - lastTime);
    totalTime += weight * d;
    return d;
  }

  double total() {
    if (on)
      return totalTime + getTime() - lastTime;
    else
      return totalTime;
  }

  double next() {
    if (!on) return 0.0;
    double t = getTime();
    double td = t - lastTime;
    totalTime += td;
    lastTime = t;
    return td;
  }

  static void reportT(double time) {
    std::cout << std::setprecision(3) << time << std::endl;
    ;
  }

  static void report_time(std::string str, double time) {
    std::cout << str << " : " << std::setprecision(3) << time << std::endl;
    ;
  }

  static void report_time_no_newline(std::string str, double time) {
    std::cout << str << " : " << std::setprecision(3) << time << ' ';
  }

  void reportTime(double time) { reportT(time); }

  void reportStop(double weight, std::string str) {
    std::cout << str << " :" << weight << ": ";
    reportTime(stop(weight));
  }

  void reportTotal() {
    double to = (totalWeight > 0.0) ? total() / totalWeight : total();
    reportTime(to);
    totalTime = 0.0;
    totalWeight = 0.0;
  }

  void reportTotal(std::string str) {
    std::cout << str << " : ";
    reportTotal();
  }

  void reportNext() { reportTime(next()); }

  void reportNext(std::string str) {
    std::cout << str << " : ";
    reportNext();
  }
};

static timer _tm;
#define timeStatement(_A, _string) \
  _tm.start();                     \
  _A;                              \
  _tm.reportNext(_string);
#define startTime() _tm.start();
#define stopTime(_weight, _str) _tm.reportStop(_weight, _str);
#define reportTime(_str) _tm.reportTotal(_str);
#define nextTime(_string) _tm.reportNext(_string);
#define nextTimeN() _tm.reportT(_tm.next());

}  // namespace elektra