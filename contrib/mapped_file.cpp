/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

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
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <iostream>
#include <string>
#include <cstring>
#include <stdexcept>

#include <cstdio>

#include "mapped_file.h"

mapped_file::mapped_file()
    : _is_open(false)
{}

mapped_file::mapped_file(const std::string& name)
    : _is_open(false)
{
    open(name);
}

mapped_file::~mapped_file()
{
    close();
}

void
mapped_file::open(const std::string& name)
{
    _is_open = map(name);
}

void
mapped_file::close(void)
{
    unmap();
}

bool
mapped_file::map(const std::string& name)
{
    int fd = ::open(name.c_str(), O_RDONLY);
    if (fd == -1)
    {
        perror("open");
        return false;
    }
    _fd = fd;

    struct stat sb;
    if (fstat(fd, &sb) == -1)
    {
        perror("fstat");
        ::close(fd);
        return false;
    }
    _length = sb.st_size;

    const char* _addr = static_cast<const char*>(
                         mmap(NULL, _length, PROT_READ, MAP_PRIVATE, fd, 0u)
                       );
    if (_addr == MAP_FAILED)
    {
        perror("mmap");
        ::close(fd);
        return false;
    }

    madvise((void*)_addr, sb.st_size, MADV_SEQUENTIAL);

    _name = name;

    _start = _addr;
    _end = _addr + _length;

    return true;

}

bool
mapped_file::unmap(void)
{
    munmap((void *)_addr, _length);
    ::close(_fd);
    return true;
}

bool
mapped_file::is_open(void) const
{
    return _is_open;
}

bool
mapped_file::end(void) const
{
    return _start && _start != _end;
}

std::string
mapped_file::get_line(void)
{
    auto old = _start;
    _start = static_cast<const char *>(memchr(_start, '\n', _end - _start));
    if (!_start)
        throw std::out_of_range("");

    _start++;

    return std::string(old, _start-1);
}

const char *
mapped_file::mem()
{
    return _start;
}



/*

int main(void)
{
    mapped_file mf("b104.mesh");

    std::cout << mf.get_line() << std::endl;
    std::cout << mf.get_line() << std::endl;
    std::cout << mf.get_line() << std::endl;
    std::cout << mf.get_line() << std::endl;

    mf.close();

    return 0;
}
*/
