#$Id: pv.glyph.py 191 2015-12-15 21:46:16Z mexas $

try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

z2_raw = ImageReader( FilePrefix='z2.raw' )

z2_raw.DataExtent = [1, 20, 1, 20, 1, 20]
z2_raw.DataByteOrder = 'LittleEndian'
z2_raw.DataScalarType = 'int'

RenderView1 = GetRenderView()
RenderView1.CameraClippingRange = [54.26289244101346, 128.9133169662742]
RenderView1.CameraFocalPoint = [9.5, 9.5, 9.5]
RenderView1.CameraParallelScale = 22.950479555412187
RenderView1.CameraPosition = [-29.88927818038211, 54.996731860164914, 72.63167352063947]
RenderView1.CameraViewUp = [-0.8592889079758231, -0.03594429642478058, -0.5102260089256964]
RenderView1.CenterOfRotation = [9.5, 9.5, 9.5]
RenderView1.ViewSize=[800,800]

DataRepresentation1 = Show()
DataRepresentation1.ScalarOpacityUnitDistance = 1.7320508075688776
DataRepresentation1.Representation = 'Outline'
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Threshold1 = Threshold()
Threshold1.Scalars = ['POINTS', 'ImageFile']
Threshold1.ThresholdRange = [1.0, 1.0]
Threshold1.AllScalars = 0

a1_ImageFile_PVLookupTable = GetLookupTableForArray( "ImageFile", 1, RGBPoints=[-60.0, 0.0, 0.0, 0.0, -35.6, 0.9019607843137255, 0.0, 0.0, -11.199999999999998, 0.9019607843137255, 0.9019607843137255, 0.0, 1.0, 1.0, 1.0, 1.0], VectorMode='Magnitude', NanColor=[0.0, 0.4980392156862745, 1.0], NumberOfTableValues=9, ColorSpace='RGB', ScalarRangeInitialized=1.0 )

a1_ImageFile_PiecewiseFunction = CreatePiecewiseFunction()

DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectMapper = 'Fixed point'
DataRepresentation2.ScalarOpacityFunction = a1_ImageFile_PiecewiseFunction
DataRepresentation2.ColorArrayName = 'ImageFile'
DataRepresentation2.ScalarOpacityUnitDistance = 2.176299194073682
DataRepresentation2.LookupTable = a1_ImageFile_PVLookupTable

Glyph1 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )
Glyph1.Vectors = ['POINTS', '']
Glyph1.GlyphTransform = "Transform2"
Glyph1.GlyphType = "Sphere"

DataRepresentation3 = Show()
DataRepresentation3.ColorArrayName = ''
DataRepresentation3.LookupTable = a1_ImageFile_PVLookupTable
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

DataRepresentation2.Visibility = 0

WriteImage('z.png')

Render()
