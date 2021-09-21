#$Id: debug.py 191 2015-12-15 21:46:16Z mexas $

from paraview.simple import *

reader=ImageReader(FilePrefix= "z9end.raw" )
reader.DataByteOrder=1
reader.DataExtent=[1,80,1,80,1,640]
reader.DataScalarType=5

view = GetRenderView()
Show()
