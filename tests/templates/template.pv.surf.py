#$Id: template.pv.surf.py 175 2015-12-15 12:31:30Z mexas $
#
# pvbatch script to view the surface of the model
#
# It will read a raw datafile specified,
# and show it as a surface with cells coloured
# by their values.
#
# Make sure to adjust the following parameters, at least.
# IMPORTANT!! PV crashes if data extents start from 1.
# So always start date extent from 0, and correspondingly
# reduce the upper extents by 1.

######################################################################72
# Adjust these parameters
ext1    = 128-1       # data extent along 1
ext2    = ext1        # data extent along 2
ext3    = ext1        # data extent along 3
infile  = "zg1.raw"   # input file name
outfile = "zg1.png"   # output file name
imsize1 = 800         # size, in pixels, of the resulting image along 1
imsize2 = 800         # size, in pixels, of the resulting image along 2
minval  = 1           # min colour value
maxval  = 3           # max colour value, equal to the number of grains

# End of adjustable parameters
######################################################################72

# define the centre of rotation (cor)
cor1 = 0.5 * ext1
cor2 = 0.5 * ext2
cor3 = 0.5 * ext3

from paraview.simple import *

reader=ImageReader(FilePrefix= infile )
reader.DataExtent=[ 0, ext1, 0, ext2, 0, ext3 ]
reader.DataByteOrder = 'LittleEndian'
reader.DataScalarType = 'int'

RenderView1 = GetRenderView()
DataRepresentation1 = Show()

dp = GetDisplayProperties(reader)
dp.LookupTable = MakeBlueToRedLT(minval,maxval)
dp.ColorAttributeType = 'POINT_DATA'
dp.ColorArrayName = 'ImageFile'
dp.Representation = "Surface"

bar = CreateScalarBar(LookupTable=dp.LookupTable, Title="grain")
bar.Position=[0.80,0.15]
GetRenderView().Representations.append(bar)

camera = GetActiveCamera()
camera.SetViewUp(-1,0,0)
camera.Azimuth(30)
camera.Elevation(30)

#camera.SetPosition(0,0,100)
#camera.Roll(-90)

RenderView1.ResetCamera()
# make white background, for papers
#RenderView1.Background = [ 1,1,1]
RenderView1.Background = [0.3176470588235294, 0.3411764705882353, 0.43137254901960786]
RenderView1.CenterAxesVisibility = 0
RenderView1.OrientationAxesVisibility = 1
RenderView1.CenterOfRotation = [ cor1, cor2, cor3 ]
RenderView1.ViewSize = [ imsize1, imsize2 ]

WriteImage( outfile )

Render()
