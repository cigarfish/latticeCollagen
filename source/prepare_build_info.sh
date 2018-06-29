#!/bin/sh
git log | head -1 | awk '{print "#define CS_BUILD_INFO_COMMIT \"",$2,"\"" }' > commit.h;
git branch | grep '*' | awk '{print "#define CS_BUILD_INFO_BRANCH \"",$2,"\"" }' >> commit.h;
diffoncommit="$(git diff | tr -d '\r' | awk '{ gsub(/\\/,"\\\\"); gsub(/"/,"\\\"") } 1' | awk -v RS="\n" -v ORS="\\\\n\"\n\"# " '1' | sed '$ d')";
cacheddiff="$(git diff --cached | tr -d '\r' | awk '{ gsub(/\\/,"\\\\"); gsub(/"/,"\\\"") } 1' | awk -v RS="\n" -v ORS="\\\\n\"\n\"# " '1' | sed '$ d')";


if [ -n "$diffoncommit" -o -n "$cacheddiff" ]; then
  /bin/echo "#define TAINTED_BY_DIFF 1" >> commit.h
  /bin/echo "const char * buildInfo_diff =" >> commit.h
  if [ -n "$cacheddiff" ]; then
      /bin/echo -n "\"# $cacheddiff" >> commit.h
  fi
  if [ -n "$diffoncommit" ]; then
     /bin/echo -n "\"# $diffoncommit" >> commit.h
  fi
  /bin/echo ";" >> commit.h
fi
