
#pragma once

#include <locale>
#include <vector>

struct csv_reader : std::ctype<char>
{
    csv_reader()
        : std::ctype<char>(get_table())
    {}

    static std::ctype_base::mask const * get_table()
    {
        static std::vector<std::ctype_base::mask>
            rc(table_size, std::ctype_base::mask());

        rc[','] = std::ctype_base::space;
        rc['\n'] = std::ctype_base::space;
        rc['\t'] = std::ctype_base::space;
        rc[' '] = std::ctype_base::space;

        return &rc[0];
    }

};


