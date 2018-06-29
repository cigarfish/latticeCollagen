/*
 * timer.h
 *
 * Andreas Buttenschoen
 * 2014
 *
 */

#include "timer.h"
#include "macros.h"
#include <chrono>

#include "spdlog/fmt/fmt.h"

namespace detail {

std::string get_unit(UNIT unit)
{
    switch (unit) {
        case UNIT::HR:
            return {" (hr)"};
        case UNIT::MIN:
            return {" (min)"};
        case UNIT::S:
            return {" (s)"};
        case UNIT::MS:
            return {" (ms)"};
        case UNIT::US:
            return {" (us)"};
        case UNIT::NS:
            return {" (ns)"};
        default:
            unreachable("Unknown time unit.");
        break;
    }
    // shut up compiler
    return {""};
}

std::string get_time_string(std::chrono::high_resolution_clock::duration duration)
{
    auto dhours = duration_cast<hours>(duration);
    auto dmin   = duration_cast<minutes>(duration-dhours);
    auto ds     = duration_cast<seconds>(duration-dmin);
    auto dms    = duration_cast<milliseconds>(duration-ds);
    auto dus    = duration_cast<microseconds>(duration-ds-dms);
    auto dns    = duration_cast<nanoseconds>(duration-ds-dms-dus);

    std::ostringstream ostr;
    UNIT unit;

    if (dhours.count()!=0) {
        unit=detail::UNIT::HR;
        ostr << fmt::format("{0:1d}", dhours.count()) << ":";
        ostr << fmt::format("{0:02d}", dmin.count()) << ":";
        ostr << fmt::format("{0:02d}", ds.count());
    } else if (dmin.count()!=0) {
        unit=detail::UNIT::MIN;
        ostr << fmt::format("{0:3d}", dmin.count()) << ":";
        ostr << fmt::format("{0:02d}", ds.count());
    } else if (ds.count()!=0) {
        unit=detail::UNIT::S;
        ostr << fmt::format("{0:3d}", ds.count()) << ".";
        ostr << fmt::format("{0:03d}", dms.count());
    } else if (dms.count()!=0) {
        unit=detail::UNIT::MS;
        ostr << fmt::format("{0:3d}", dms.count()) << ".";
        ostr << fmt::format("{0:03d}", dus.count());
    } else if (dus.count()!=0) {
        unit=detail::UNIT::US;
        ostr << fmt::format("{0:3d}", dus.count()) << ".";
        ostr << fmt::format("{0:03d}", dns.count());
    } else {
        unit=detail::UNIT::NS;
        ostr << dns.count();
    }

    ostr << get_unit(unit);

    return ostr.str();
}

} // end namespace
