#$Id: template.pv.cracks.py 175 2015-12-15 12:31:30Z mexas $
#
# pvbatch script to show threshold values, in this
# case cracks, i.e. cells with states between some
# lower and upper values.
#
# The script reads the raw binary file specified.
# The cracks are shown as surfaces with a uniform colour.
#
# Make sure to adjust the following parameters, at least.
# IMPORTANT!! PV crashes if data extents start from 1.
# So always start date extent from 0, and correspondingly
# reduce the upper extents by 1.

# CGPACK cell states
#                  cgca_gb_state_intact: 2
#               cgca_gb_state_fractured: 1
#                     cgca_intact_state: 0
#             cgca_clvg_state_100_flank: -1
#              cgca_clvg_state_100_edge: -2
#             cgca_clvg_state_110_flank: -3
#              cgca_clvg_state_110_edge: -4
#             cgca_clvg_state_111_flank: -5
#              cgca_clvg_state_111_edge: -6
# end CGPACK cell states

######################################################################72
# Adjust these parameters
#1000000 0.9991665 464 (10.00000,10.00000,9.989224)
# dataset sizes for ParaView 2*4640,  4635
#Total cells in the model  9.98E+10
ext1    = 4640       # data extent along 1
ext2    = 4640       # data extent along 2
ext3    = 4635       # data extent along 3
ffile   = "./large/250/zf3.raw" # fracture file
imsize1 = 1200         # size, in pixels, of the resulting image along 1
imsize2 = 1200         # size, in pixels, of the resulting image along 2
# End of adjustable parameters
######################################################################72

# define the centre of rotation (cor)
cor1 = 0.5 * ext1
cor2 = 0.5 * ext2
cor3 = 0.5 * ext3

from paraview.simple import *

# the extents start from zero, so need to lower
# the upper extents by 1
cracks = ImageReader( FilePrefix= ffile )
cracks.DataExtent=[ 0, ext1-1, 0, ext2-1, 0, ext3-1 ]
cracks.DataByteOrder = 'LittleEndian'
cracks.DataScalarType = 'int'

RenderView1 = GetRenderView()
DataRepresentation1 = Show()

DataRepresentation1.Representation = 'Outline'
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5]

# grain boundaries, cell state 2

#Threshold1 = Threshold()
#Threshold1.Scalars = ['POINTS', 'ImageFile']
#Threshold1.ThresholdRange = [ 2, 2 ]
#Threshold1.AllScalars = 0

#DataRepresentation2 = Show()
#DataRepresentation2.ScalarOpacityUnitDistance = 1.0
#DataRepresentation2.SelectionPointFieldDataArrayName = 'ImageFile'
#DataRepresentation2.DiffuseColor = [ 1, 0, 1 ]

camera = GetActiveCamera()
camera.SetViewUp(-1,0,0)
camera.Azimuth(30)
camera.Elevation(30)

RenderView1.ResetCamera()

# gradient background colour
RenderView1.UseGradientBackground = 1
RenderView1.Background2 = [0.0, 0.0, 0.16470588235294117]
RenderView1.Background = [0.3215686274509804, 0.3411764705882353, 0.43137254901960786]

RenderView1.CenterAxesVisibility = 0
RenderView1.OrientationAxesVisibility = 1
RenderView1.CenterOfRotation = [ cor1, cor2, cor3 ]
RenderView1.CameraFocalPoint = [ cor1, cor2, cor3 ]
RenderView1.ViewSize = [ imsize1, imsize2 ]

# this is effectively a zoom. A positive value
# will zoom in, a negative will zoom out
# The value of 1 will put the whole model
# into view.
#RenderView1.CameraViewAngle = 10
RenderView1.CameraViewAngle = 30

#DataRepresentation2.ScalarOpacityFunction = a1_ImageFile_PiecewiseFunction
#DataRepresentation2.LookupTable = a1_ImageFile_PVLookupTable

# do all crack states from the main dataset
SetActiveSource( cracks )

# (100) cracks
cracks100 = Threshold()
cracks100.Scalars = ['POINTS', 'ImageFile']
cracks100.ThresholdRange = [ -2, -1 ]
cracks100.AllScalars = 0

DataRepresentation3 = Show()
DataRepresentation3.ScalarOpacityUnitDistance = 1.0
DataRepresentation3.SelectionPointFieldDataArrayName = 'ImageFile'
#DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]
#DataRepresentation3.ScaleFactor = -2.0000000000000002e+298
DataRepresentation3.DiffuseColor = [1.0, 1.0, 0.0]

# (110) cracks
SetActiveSource( cracks )
cracks110 = Threshold()
cracks110.Scalars = ['POINTS', 'ImageFile']
cracks110.ThresholdRange = [ -4, -3 ]
cracks110.AllScalars = 0

DataRepresentation4 = Show()
DataRepresentation4.ScalarOpacityUnitDistance = 1.0
DataRepresentation4.SelectionPointFieldDataArrayName = 'ImageFile'
#DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]
#DataRepresentation4.ScaleFactor = -2.0000000000000002e+298
DataRepresentation4.DiffuseColor = [0.0, 1.0, 0.5]

#DataRepresentation3.ColorArrayName = 'ImageFile'
#DataRepresentation3.SelectMapper = 'Fixed point'

# 1 is to show, 0 not to show
# data2 is GB
# data3 is cracks
# data4 is grains microstructure

#DataRepresentation2.Opacity = 0.1
#WriteImage( "crgb.png" )

#DataRepresentation2.Opacity = 1
#DataRepresentation3.Visibility = 0
#WriteImage( "gb.png" )

#DataRepresentation2.Visibility = 0
#DataRepresentation3.Visibility = 1
WriteImage( "cr.png" )

RenderView1.ResetCamera()

Render()
