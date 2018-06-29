#include "StringUtils.h"

#include <sstream>
#include <algorithm>

std::vector<std::string> splitString(std::string str, char sep)
{
    std::vector<std::string> vecString;
    std::string item;

    std::stringstream stringStream(str);

    while (std::getline(stringStream, item, sep))
        vecString.push_back(item);

    return vecString;
}

std::string& lower_case(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

std::string lower_case(const std::string& str)
{
    std::string ret;
    std::transform(str.begin(), str.end(), ret.begin(), ::tolower);
    return ret;
}