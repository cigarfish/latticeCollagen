#!/bin/bash -x

PROJECTS="CSCore CSGUI CellSysUI"

for proj in $PROJECTS; do
  pushd $proj.xcodeproj
  if [ -e project.pbxproj.* ]; then
      sed 's#path = "\(.*\)";#path = "\1"; sourceTree = SOURCE_ROOT;#g' project.pbxproj.* > project.pbxproj
  fi
  popd
done

exit 0;
