#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <vector>

#define STRING_REMOVE_CHAR(str, ch) str.erase(std::remove(str.begin(), str.end(), ch), str.end())

std::vector<std::string> splitString(std::string str, char sep = ',');

std::string& lower_case(std::string& str);

std::string lower_case(const std::string& str);

#endif