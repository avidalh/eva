#!/bin/bash
echo Batch Analisys Starts...
date

#find /Volumes/UNTITLED\ 1/grabaciones/tenerife/smr/11/1*/* -name *.gps | xargs -I{} ./AstTool_0300.py -f{} -d1 -b
#find ../../Asterix/SMR/06/16/* -name *.gps | xargs -I{} ./AstTool_0400.py -f{} -d1 -b

find /Volumes/UNTITLED\ 1/grabaciones/tenerife/smr/2014/12/*/* -name *.gps | xargs -I{} ./eva_0500.py -f{} -l1 -b -d -m -i 

echo Batch Analisys Ends...
date

