#$Id: gv-bug.awk 170 2015-11-04 11:12:30Z mexas $

# There seem to be a bug in graphviz fig filter.
# GV writes a rectangle (polygon) the size of the
# whole image as the first object.
# This object then leads to the creation of the
# bounding line around the complete image.
# This script removes the first polygon.

BEGIN{done=0;started=0}
{ if ( started==1 && $0~/^ / ) { $0="#***GV BUG?*** " $0; done=1 } else {started=0}
  if ( done==0 && $0~/^2 3 / ) { $0="#***GV BUG?*** " $0; started=1 }
print }
