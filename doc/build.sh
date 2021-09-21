#!/bin/sh
#$Id: build.sh 525 2018-03-19 21:54:26Z mexas $
#  AUTHOR
#    Anton Shterenlikht
#  COPYRIGHT
#    See LICENSE
#  DESCRIPTION
#    This shell scripts builds all ROBOdoc documentation and saves it
#    to a directory. This directory can then be moved/copied to the
#    main casup doc directory for upload to sf.net.
#
#    Mare sure to run the script from casup/head/doc directory.
#    All paths are relative to that location.

STATS=stats
ROBO_DIR=robodoc
CASUP=casup

# Move to casup/head

cd ..

# Create file $STATS with code statistics

echo "!*robodoc*e* CASUP/stats" > $STATS
echo "!  NAME"    >> $STATS
echo "!    stats" >> $STATS
echo "!  SOURCE"  >> $STATS

cloc --by-file *f90 >> $STATS

echo "!*roboend*" >> $STATS

# Build all robodoc files under casup/head. This is a recursive
# build down the directory tree, .i.e. casup/head/tests will
# be processed too.

# First make a new directory $ROBO_DIR to store all html files.
# If the directlry already exists, delete it first.

if [ -d $ROBO_DIR ]; then
	rm -rf $ROBO_DIR
fi
mkdir $ROBO_DIR

robodoc --html --multidoc --doc $ROBO_DIR/
chmod 755 $ROBO_DIR $ROBO_DIR/tests
robodoc --ascii --singledoc --documenttitle "CASUP source code and tests" --doc $ROBO_DIR/$CASUP
robodoc --latex --altlatex --singledoc --documenttitle "CASUP source code and tests" --doc $ROBO_DIR/$CASUP

# Remove STATS now
rm $STATS

# Build PDF from latex src
cd $ROBO_DIR
pdflatex $CASUP
pdflatex $CASUP
pdflatex $CASUP
rm $CASUP.aux $CASUP.idx $CASUP.log $CASUP.tex $CASUP.toc
