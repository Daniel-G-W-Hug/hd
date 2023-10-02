#pragma once

//
// provide a stop_watch to measure execution time
//

#include "date/date.h"
#include <algorithm> // std::min
#include <chrono>
#include <sstream>
#include <string>
#include <vector>

namespace hd {

std::string now_as_str()
{
    // use Howard Hinnants data.h to print current time
    auto now = std::chrono::system_clock::now();
    auto today = date::floor<date::days>(now);

    using date::operator<<; // for stream output of today
    std::stringstream ss;
    ss << today << ' ' << date::make_time(now - today) << " UTC";
    return ss.str();

    // #include "date/tz.h"
    // #include <iostream>
    // std::cout << date::zoned_time{date::current_zone(),
    // std::chrono::system_clock::now()};
}

enum class time_in { seconds,
                     milliseconds,
                     microseconds,
                     nanoseconds };

// Usage of class stop_watch:
//
// Call start()-stop()-pairs at least once.
// Corresponding time differences for each pair will be accumulated in elapsed_time().
//
// split() will effectively call stop() and start() using the same time point, i.e. it
// should occur between calls to start() and stop() when the split time is required to use
// for output with elapsed_time().
//
// elapsed_time() returns the accumulated time since the last call to start(). Can be
// called after stop() or split() (or otherwise returns 0).
//
// After using the stop watch a call to reset() does what is promises. Without calling
// reset(), further calls to elapsed_time() will continue to accumulate all previous and
// all new calls to start() or stop() until that point time.

class stop_watch {

    std::vector<std::chrono::steady_clock::time_point> start_time{};
    std::vector<std::chrono::steady_clock::time_point> end_time{};
    int start_cnt{0};
    int stop_cnt{0};

  public:
    void start(); // start time
    void split(); // split time (interim time between start and stop)
    void stop();  // stop time
    int elapsed_time(time_in t_in);
    void reset();
};

void stop_watch::start()
{
    start_time.push_back(std::chrono::steady_clock::now());
    ++start_cnt;
}

void stop_watch::split()
{
    auto now = std::chrono::steady_clock::now();
    // end of current interval (initiated by start() or split())
    end_time.push_back(now);
    ++stop_cnt;
    // begin of new interval (to be finished with stop() or split())
    start_time.push_back(now);
    ++start_cnt;
}

void stop_watch::stop()
{
    end_time.push_back(std::chrono::steady_clock::now());
    ++stop_cnt;
}

int stop_watch::elapsed_time(time_in t_in)
{

    using namespace std::chrono;

    int complete_measurements = std::min(start_cnt, stop_cnt);

    if (complete_measurements >= 1) { // return values for complete measurements

        int time_difference = 0;

        switch (t_in) {
        case time_in::seconds:
            for (int i = 0; i < complete_measurements; ++i) {
                time_difference +=
                    duration_cast<seconds>(end_time[i] - start_time[i]).count();
            }
            return time_difference;
            break;
        case time_in::milliseconds:
            for (int i = 0; i < complete_measurements; ++i) {
                time_difference +=
                    duration_cast<milliseconds>(end_time[i] - start_time[i]).count();
            }
            return time_difference;
            break;
        case time_in::microseconds:
            for (int i = 0; i < complete_measurements; ++i) {
                time_difference +=
                    duration_cast<microseconds>(end_time[i] - start_time[i]).count();
            }
            return time_difference;
            break;
        case time_in::nanoseconds:
            for (int i = 0; i < complete_measurements; ++i) {
                time_difference +=
                    duration_cast<nanoseconds>(end_time[i] - start_time[i]).count();
            }
            return time_difference;
            break;
        }
    }
    else
        return 0; // otherwise return 0
}

void stop_watch::reset()
{
    start_time.clear();
    end_time.clear();
    start_cnt = 0;
    stop_cnt = 0;
}

} // namespace hd