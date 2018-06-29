#ifndef FILES_H
#define FILES_H

#include <sys/stat.h>
#include <string>

inline bool FileExists(const std::string& filename)
{
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

#endif
