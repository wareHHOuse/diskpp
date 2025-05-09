//
//  display_settings.cpp
//  acoustics
//
//  Created by Omar Durán on 4/8/20.
//

#include <unistd.h>

std::ostream& bold(std::ostream& os) { os << "\x1b[1m"; return os; }
std::ostream& nobold(std::ostream& os) { os << "\x1b[21m"; return os; }

std::ostream& underline(std::ostream& os) { os << "\x1b[4m"; return os; }
std::ostream& nounderline(std::ostream& os) { os << "\x1b[24m"; return os; }

std::ostream& blink(std::ostream& os) { os << "\x1b[5m"; return os; }
std::ostream& noblink(std::ostream& os) { os << "\x1b[25m"; return os; }

std::ostream& reset(std::ostream& os) { os << "\x1b[0m"; return os; }
std::ostream& erase_line(std::ostream& os) { os << "\x1b[0K"; return os; }

std::ostream& red(std::ostream& os) { os << "\x1b[31m"; return os; }
std::ostream& green(std::ostream& os) { os << "\x1b[32m"; return os; }
std::ostream& yellow(std::ostream& os) { os << "\x1b[33m"; return os; }
std::ostream& blue(std::ostream& os) { os << "\x1b[34m"; return os; }
std::ostream& magenta(std::ostream& os) { os << "\x1b[35m"; return os; }
std::ostream& cyan(std::ostream& os) { os << "\x1b[36m"; return os; }
std::ostream& nocolor(std::ostream& os) { os << "\x1b[39m"; return os; }
