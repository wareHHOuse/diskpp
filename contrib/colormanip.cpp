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

#include <time.h>
#include "colormanip.h"

/* COLORS */
std::ostream& red(std::ostream& os) { os << "\x1b[31m"; return os; }
std::ostream& green(std::ostream& os) { os << "\x1b[32m"; return os; }
std::ostream& yellow(std::ostream& os) { os << "\x1b[33m"; return os; }
std::ostream& blue(std::ostream& os) { os << "\x1b[34m"; return os; }
std::ostream& magenta(std::ostream& os) { os << "\x1b[35m"; return os; }
std::ostream& cyan(std::ostream& os) { os << "\x1b[36m"; return os; }
std::ostream& nocolor(std::ostream& os) { os << "\x1b[39m"; return os; }

/* BACKGROUND COLORS */
std::ostream& bgred(std::ostream& os) { os << "\x1b[41m"; return os; }
std::ostream& bggreen(std::ostream& os) { os << "\x1b[42m"; return os; }
std::ostream& bgyellow(std::ostream& os) { os << "\x1b[43m"; return os; }
std::ostream& bgblue(std::ostream& os) { os << "\x1b[44m"; return os; }
std::ostream& bgmagenta(std::ostream& os) { os << "\x1b[45m"; return os; }
std::ostream& bgcyan(std::ostream& os) { os << "\x1b[46m"; return os; }
std::ostream& nobg(std::ostream& os) { os << "\x1b[49m"; return os; }

/* BOLD (boldoff widely unsupported!) */
std::ostream& bold(std::ostream& os) { os << "\x1b[1m"; return os; }
std::ostream& nobold(std::ostream& os) { os << "\x1b[21m"; return os; }

/* UNDERLINE */
std::ostream& underline(std::ostream& os) { os << "\x1b[4m"; return os; }
std::ostream& nounderline(std::ostream& os) { os << "\x1b[24m"; return os; }

/* BLINK */
std::ostream& blink(std::ostream& os) { os << "\x1b[5m"; return os; }
std::ostream& noblink(std::ostream& os) { os << "\x1b[25m"; return os; }

/* RESET */
std::ostream& reset(std::ostream& os) { os << "\x1b[0m"; return os; }


/* TIME */
std::ostream& time_now(std::ostream& os)
{
    time_t      rawtime;
    struct tm   *timeinfo;
    char        buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime (buffer,80,"[%D %T] ",timeinfo);

    os << buffer;
    return os;
}
