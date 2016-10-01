/*
 * This source file is part of EMT, the ElectroMagneticTool.
 *
 * Copyright (C) 2013-2015, Matteo Cicuttin - matteo.cicuttin@uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <map>

#include <sys/resource.h>

using namespace std::chrono;

class timecounter
{
    steady_clock::time_point            m_start, m_stop;
    bool                                m_is_running;

public:
                                        timecounter();

    steady_clock::time_point            tic();
    steady_clock::time_point            toc();
    duration<double>                    elapsed() const;
    bool                                is_running() const;
    double                              to_double() const;
};


class timecounter_new
{
    struct rusage m_start, m_stop;

public:
    timecounter_new()
    {}

    void tic()
    {
        getrusage(RUSAGE_SELF, &m_start);
    }

    void toc()
    {
        getrusage(RUSAGE_SELF, &m_stop);
    }

    double get_usertime() const
    {
        double start, stop;
        start = m_start.ru_utime.tv_sec + double(m_start.ru_utime.tv_usec)/1e6;
        stop = m_stop.ru_utime.tv_sec + double(m_stop.ru_utime.tv_usec)/1e6;
        return stop - start;
    }

    double get_systime() const
    {
        double start, stop;
        start = m_start.ru_stime.tv_sec + double(m_start.ru_stime.tv_usec)/1e6;
        stop = m_stop.ru_stime.tv_sec + double(m_stop.ru_stime.tv_usec)/1e6;
        return stop - start;
    }

    double to_double() const
    {
        return get_systime() + get_usertime();
    }
};

std::ostream&
operator<<(std::ostream& os, const timecounter_new& tc)
{
    os << tc.to_double();

    return os;
}


class time_profiler
{
    typedef typename steady_clock::time_point                   time_point;

    std::map<std::string, std::pair<time_point, time_point>>    m_timings;
    std::map<std::string, size_t>                               m_calls;
    std::map<std::string, duration<double>>                     m_durations;

public:
            time_profiler();

    void    enter(const std::string&);
    void    leave(const std::string&);
    void    report(void);
};

#ifdef PROFILING
    #define PROFILER_ENTER(profiler, counter) { profiler.tic(counter); }
    #define PROFILER_LEAVE(profiler, counter) { profiler.toc(counter); }
#else
    #define PROFILER_ENTER(profiler, counter)
    #define PROFILER_LEAVE(profiler, counter)
#endif







std::ostream& operator<<(std::ostream&, const timecounter&);
