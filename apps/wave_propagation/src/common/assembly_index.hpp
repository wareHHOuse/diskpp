//
//  assembly_index.hpp
//  acoustics
//
//  Created by Omar Durán on 5/6/20.
//

#pragma once
#ifndef assembly_index_hpp
#define assembly_index_hpp

#include <iostream>

class assembly_index
{
    size_t  idx;
    bool    assem;

public:
    assembly_index(size_t i, bool as)
        : idx(i), assem(as)
    {}

    operator size_t() const
    {
        if (!assem)
            throw std::logic_error("Invalid assembly_index");

        return idx;
    }

    bool assemble() const
    {
        return assem;
    }
    
    size_t vidx() const
    {
        return idx;
    }

    friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
    {
        os << "(" << as.idx << "," << as.assem << ")";
        return os;
    }
};

#endif /* assembly_index_hpp */
