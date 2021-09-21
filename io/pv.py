# A basic script for pvbatch
# It will read a raw datafile specified,
# and show it as a surface with cells coloured
# by their values.
# Make sure to adjust the following parameters, at least.

##########################################################
# Adjust these parameters
ext1=1024 - 1
ext2=ext1
ext3=ext1
minval=1
maxval=512
infile="native.dat"
outfile="z.png"
# End of adjustable parameters
##########################################################

from paraview.simple import *

reader=ImageReader(FilePrefix= infile )
reader.DataByteOrder=1
reader.DataExtent=[ 0, ext1, 0, ext2, 0, ext3 ]
reader.DataScalarType=6

view = GetActiveView()
if not view:
 view = CreateRenderView()
view.ViewSize=[800,800]
view.Background=[0.3249412, 0.34902, 0.427451]

Show()

dp = GetDisplayProperties(reader)
dp.LookupTable = MakeBlueToRedLT(minval,maxval)
dp.ColorAttributeType = 'POINT_DATA'
dp.ColorArrayName = 'ImageFile'
dp.Representation = "Surface"

bar = CreateScalarBar(LookupTable=dp.LookupTable, Title="grain")
#bar.Position=[0.80,0.15]

GetRenderView().Representations.append(bar)

camera = GetActiveCamera()
camera.SetViewUp(-1,0,0)
camera.Azimuth(30)
camera.Elevation(30)

#camera.SetPosition(0,0,100)
#camera.Roll(-90)

ResetCamera()

WriteImage( outfile )

Render()
