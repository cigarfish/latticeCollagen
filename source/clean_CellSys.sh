#!/bin/sh -x

rm -rf CSCore.xcodeproj CSGUI.xcodeproj source.xcodeproj CellSysUI.xcodeproj moc_*

qmake -recursive -spec macx-xcode

sh correct-xcodeproj-gg.sh

open source.xcodeproj
