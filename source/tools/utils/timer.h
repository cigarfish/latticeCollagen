/*
 * timer.h
 *
 * Andreas Buttenschoen
 * 2014
 *
 */

#ifndef TIMER_H
#define TIMER_H

#define TIMER_H_VERSION 0.2

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

#include "macros.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::hours;
using std::chrono::minutes;
using std::chrono::seconds;
using std::chrono::milliseconds;
using std::chrono::microseconds;
using std::chrono::nanoseconds;

namespace detail {

enum class UNIT : unsigned short { HR, MIN, S, MS, US, NS };

std::string get_unit(UNIT unit);
std::string get_time_string(std::chrono::high_resolution_clock::duration duration);

} // end namespace

class Chronos {
public:
    Chronos()
        : start(high_resolution_clock::now()), times()
    {}

    Chronos(const Chronos& other)
        : start(other.start), times(other.times)
    {}

    Chronos& operator=(Chronos other)
    {
        std::swap(start, other.start);
        std::swap(times, other.times);
        return *this;
    }

    ~Chronos()
    {}

    void lap()
    {
        auto current_time = high_resolution_clock::now();
        times.push_back(current_time);
    }

    std::chrono::high_resolution_clock::duration time() const
    {
        return high_resolution_clock::now() - start;
    }

    std::chrono::high_resolution_clock::duration compute_recent_lap() const
    {
        std::size_t sz = times.size();
        return times[sz-1] - times[sz-2];
    }

    std::string recent_lap() const
    {
        return detail::get_time_string(compute_recent_lap());
    }

    std::string recent_lap(const unsigned int numberOfSubIntervals)
    {
        if (numberOfSubIntervals > 0)
            return detail::get_time_string(compute_recent_lap() / numberOfSubIntervals);
        else
            return detail::get_time_string(time());
    }

    std::string getTime() const
    {
        return detail::get_time_string(time());
    }

    void reset()
    {
        start = high_resolution_clock::now();
        times.clear();
        times.push_back(start);
    }

    void print(std::ostream& os)
    {
        compute(os);
    }

private:

    void compute(std::ostream& os)
    {
        auto end = high_resolution_clock::now();

        auto dhours = duration_cast<hours>(end-start);
        auto dmin = duration_cast<minutes>(end-start-dhours);
        auto ds = duration_cast<seconds>(end-start-dmin);
        auto dms = duration_cast<milliseconds>(end-start-ds);
        auto dus = duration_cast<microseconds>(end-start-ds-dms);
        auto dns = duration_cast<nanoseconds>(end-start-ds-dms-dus);

        std::ostringstream ostr;
        detail::UNIT unit;

        ostr << "Program Execution Time: ";
        if (dhours.count()!=0) {
            unit=detail::UNIT::HR;
            ostr << dhours.count() << ":";
            ostr << dmin.count() << ":";
            ostr << ds.count();
        } else if (dmin.count()!=0) {
            unit=detail::UNIT::MIN;
            ostr << dmin.count() << ":";
            ostr << ds.count();
        } else if (ds.count()!=0) {
            unit=detail::UNIT::S;
            ostr << ds.count() << ".";
            ostr << dms.count();
        } else if (dms.count()!=0) {
            unit=detail::UNIT::MS;
            ostr << dms.count() << ".";
            ostr << dus.count();
        } else if (dus.count()!=0) {
            unit=detail::UNIT::US;
            ostr << dus.count() << ".";
            ostr << dns.count();
        } else {
            unit=detail::UNIT::NS;
            ostr << dns.count();
        }

        ostr << detail::get_unit(unit);
        ostr << std::endl;

        os << ostr.str();
    }

    high_resolution_clock::time_point start;
    std::vector<high_resolution_clock::time_point> times;
};

#endif
