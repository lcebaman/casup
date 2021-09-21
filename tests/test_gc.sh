#!/bin/sh
#
#$Id: test_gc.sh 9 2014-12-01 09:55:21Z mexas $
#
# The only argument to this program is the common
# part of the grain connectivity files, created
# by tests AAY, AAZ or ABA.
cat $1* | sort -n -k1 -k2 | uniq > tmp
./test_gc.xnonca < tmp
rm tmp
