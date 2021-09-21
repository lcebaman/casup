#$Id: pv.4files.py 191 2015-12-15 21:46:16Z mexas $

try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

#parameters
tmin=1
tmax=1
ext1=40
ext2=ext1
ext3=ext1
imgsize=1000

# input and output files
in1='z0.raw'
in2='z10.raw'
in3='z20.raw'
in4='z30.raw'

out1="z1.png"
out2="z2.png"
out3="z3.png"
out4="z4.png"

#**********************************************************************
# first file
#**********************************************************************

z1_raw = ImageReader( FilePrefix=in1 )
z1_raw.DataExtent = [1,ext1,1,ext2,1,ext3]
z1_raw.DataByteOrder = 'LittleEndian'
z1_raw.DataScalarType = 'int'

RenderView1 = GetRenderView()
RenderView1.CameraClippingRange = [116.40010427154998, 258.29444178717074]
RenderView1.CameraFocalPoint = [19.5, 19.5, 19.5]
RenderView1.CameraParallelScale = 47.108879087425024
RenderView1.CameraPosition = [-39.95937860419534, 85.12019188293596, 175.09382644680005]
RenderView1.CameraViewUp = [-0.927564936038639, -0.29417702924143097, -0.23039784048103]
RenderView1.CenterOfRotation = [19.5, 19.5, 19.5]
RenderView1.ViewSize=[imgsize, imgsize]

DataRepresentation1 = Show()
DataRepresentation1.ScalarOpacityUnitDistance = 1.7320508075688776
DataRepresentation1.Representation = 'Outline'
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Threshold1 = Threshold()
Threshold1.Scalars = ['POINTS', 'ImageFile']
Threshold1.ThresholdRange = [tmin,tmax]
Threshold1.AllScalars = 0

a1_ImageFile_PVLookupTable = GetLookupTableForArray( "ImageFile", 1 )

DataRepresentation2 = Show()
DataRepresentation2.ColorArrayName = 'ImageFile'
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.LookupTable = a1_ImageFile_PVLookupTable
DataRepresentation2.Representation = 'Surface With Edges'
DataRepresentation2.SelectMapper = 'Fixed point'
DataRepresentation2.ScalarOpacityFunction = []
DataRepresentation2.ScalarOpacityUnitDistance = 1.7320508075688776

WriteImage( out1 )

#**********************************************************************
# second file
#**********************************************************************

z2_raw = ImageReader( FilePrefix=in2 )
z2_raw.DataExtent = [1,ext1,1,ext2,1,ext3]
z2_raw.DataByteOrder = 'LittleEndian'
z2_raw.DataScalarType = 'int'

DataRepresentation3 = Show()
DataRepresentation3.ScalarOpacityUnitDistance = 1.7320508075688776
DataRepresentation3.Representation = 'Outline'
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Delete(z1_raw)

Threshold2 = Threshold()
Threshold2.Scalars = ['POINTS', 'ImageFile']
Threshold2.ThresholdRange = [tmin,tmax]
Threshold2.AllScalars = 0

DataRepresentation4 = Show()
DataRepresentation4.ColorArrayName = ''
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation4.LookupTable = a1_ImageFile_PVLookupTable
DataRepresentation4.Representation = 'Surface With Edges'
DataRepresentation4.SelectMapper = 'Fixed point'
#DataRepresentation4.ScalarOpacityFunction = a1_ImageFile_PiecewiseFunction
DataRepresentation4.ScalarOpacityFunction = []
DataRepresentation4.ScalarOpacityUnitDistance = 2.0537331866470296

a1_ImageFile_PVLookupTable.RGBPoints = [-55.0, 0.0, 0.0, 0.0, -32.599999999999994, 0.9019607843137255, 0.0, 0.0, -10.199999999999996, 0.9019607843137255, 0.9019607843137255, 0.0, 1.0, 1.0, 1.0, 1.0]

WriteImage( out2 )

#**********************************************************************
# third file
#**********************************************************************

z3_raw = ImageReader( FilePrefix=in3 )
z3_raw.DataExtent = [1,ext1,1,ext2,1,ext3]
z3_raw.DataByteOrder = 'LittleEndian'
z3_raw.DataScalarType = 'int'

DataRepresentation5 = Show()
DataRepresentation5.ScalarOpacityUnitDistance = 1.7320508075688776
DataRepresentation5.Representation = 'Outline'
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Delete(z2_raw)

Threshold3 = Threshold()
Threshold3.Scalars = ['POINTS', 'ImageFile']
Threshold3.ThresholdRange = [tmin,tmax]
Threshold3.AllScalars = 0

DataRepresentation6 = Show()
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation6.SelectMapper = 'Fixed point'
DataRepresentation6.ScalarOpacityFunction = []
DataRepresentation6.ColorArrayName = 'ImageFile'
DataRepresentation6.ScalarOpacityUnitDistance = 2.2335530257283205
DataRepresentation6.LookupTable = a1_ImageFile_PVLookupTable
DataRepresentation6.Representation = 'Surface With Edges'
DataRepresentation6.ColorArrayName = ''

a1_ImageFile_PVLookupTable.RGBPoints = [-60.0, 0.0, 0.0, 0.0, -35.599999999999994, 0.9019607843137255, 0.0, 0.0, -11.2, 0.9019607843137255, 0.9019607843137255, 0.0, 1.0, 1.0, 1.0, 1.0]

WriteImage( out3 )

#**********************************************************************
# last file
#**********************************************************************

z9end_raw = ImageReader( FilePrefix=in4 )
z9end_raw.DataExtent = [1,ext1,1,ext2,1,ext3]
z9end_raw.DataByteOrder = 'LittleEndian'
z9end_raw.DataScalarType = 'int'

DataRepresentation7 = Show()
DataRepresentation7.ScalarOpacityUnitDistance = 1.7320508075688776
DataRepresentation7.Representation = 'Outline'
DataRepresentation7.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Delete(z3_raw)

DataRepresentation7.Representation = 'Surface With Edges'

WriteImage( out4 )

Render()
