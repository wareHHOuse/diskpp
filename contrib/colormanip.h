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

#include <iostream>

std::ostream& red(std::ostream& os);
std::ostream& green(std::ostream& os);
std::ostream& yellow(std::ostream& os);
std::ostream& blue(std::ostream& os);
std::ostream& magenta(std::ostream& os);
std::ostream& cyan(std::ostream& os);
std::ostream& nocolor(std::ostream& os);

std::ostream& bgred(std::ostream& os);
std::ostream& bggreen(std::ostream& os);
std::ostream& bgyellow(std::ostream& os);
std::ostream& bgblue(std::ostream& os);
std::ostream& bgmagenta(std::ostream& os);
std::ostream& bgcyan(std::ostream& os);
std::ostream& nobg(std::ostream& os);

std::ostream& bold(std::ostream& os);
std::ostream& nobold(std::ostream& os);
std::ostream& underline(std::ostream& os);
std::ostream& nounderline(std::ostream& os);
std::ostream& blink(std::ostream& os);
std::ostream& noblink(std::ostream& os);

std::ostream& reset(std::ostream& os);

std::ostream& time_now(std::ostream& os);
