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

#include "timecounter.h"

using namespace std::chrono;

timecounter::timecounter()
    : m_is_running(false)
{}

steady_clock::time_point
timecounter::tic()
{
    m_start = steady_clock::now();
    m_is_running = true;
    return m_start;
}

steady_clock::time_point
timecounter::toc()
{
    m_stop = steady_clock::now();
    m_is_running = false;
    return m_stop;
}

duration<double>
timecounter::elapsed() const
{
    if (m_is_running)
    {
        auto now = steady_clock::now();
        return now - m_start;
    }

    return m_stop - m_start;
}

bool
timecounter::is_running() const
{
    return m_is_running;
}

double
timecounter::to_double() const
{
    duration<double> time_span = duration_cast<duration<double>>(elapsed());
    return time_span.count();
}

std::ostream&
operator<<(std::ostream& os, const timecounter& tc)
{
    duration<double> time_span = duration_cast<duration<double>>(tc.elapsed());
    if (tc.is_running())
        os << time_span.count() << " and counting";
    else
        os << time_span.count();

    return os;
}
