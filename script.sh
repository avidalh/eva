#!/bin/bash
echo Batch Analisys Starts...
date

#find /Volumes/UNTITLED\ 1/grabaciones/tenerife/smr/11/1*/* -name *.gps | xargs -I{} ./AstTool_0300.py -f{} -d1 -b
#find ../../Asterix/SMR/06/16/* -name *.gps | xargs -I{} ./AstTool_0400.py -f{} -d1 -b

find /Volumes/UNTITLED\ 1/grabaciones/tenerife/smr/2015/02/*/* -name *.gps | xargs -I{} ./eva.py -f{} -l1 -b -d -m -i 

echo Batch Analisys Ends...
date


