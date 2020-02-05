# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-2 replay file
# Internal Version: 2013_07_18-13.24.06 126428
# Run by dp44dafy on Wed Feb  5 19:44:46 2020
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=317.129180908203, 
    height=211.592514038086)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
Mdb()
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(-5.0, 5.0), point2=(-7.5, 11.25))
s.undo()
s.rectangle(point1=(-5.0, -5.0), point2=(5.0, 5.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=180.956, 
    farPlane=196.168, width=73.264, height=36.5404, cameraPosition=(-0.366163, 
    -2.62832, 188.562), cameraTarget=(-0.366163, -2.62832, 0))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.4196, 
    farPlane=30.1489, width=23.6071, height=11.7741, viewOffsetX=1.03918, 
    viewOffsetY=0.274789)
del mdb.models['Model-1'].parts['Part-1']
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.rectangle(point1=(-5.0, -5.0), point2=(5.0, 5.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=187.642, 
    farPlane=189.481, width=8.85802, height=4.41794, cameraPosition=(-4.20775, 
    -4.75796, 188.562), cameraTarget=(-4.20775, -4.75796, 0))
s1.CircleByCenterPerimeter(center=(-5.0, -5.0), point1=(0.333, 0.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=186.074, 
    farPlane=191.05, width=33.1711, height=16.5441, cameraPosition=(-3.5008, 
    -4.5971, 188.562), cameraTarget=(-3.5008, -4.5971, 0))
s1.autoTrimCurve(curve1=g[6], point1=(-11.8669528961182, -2.69607591629028))
s1.autoTrimCurve(curve1=g[2], point1=(-5.09435510635376, -2.34669828414917))
s1.autoTrimCurve(curve1=g[5], point1=(-1.8233790397644, -4.87454986572266))
session.viewports['Viewport: 1'].view.setValues(nearPlane=185.118, 
    farPlane=192.006, width=45.9115, height=22.8984, cameraPosition=(-3.39176, 
    -4.8737, 188.562), cameraTarget=(-3.39176, -4.8737, 0))
s1.undo()
s1.undo()
s1.undo()
s1.undo()
s1.CircleByCenterPerimeter(center=(-5.0, -5.0), point1=(0.333, 0.0))
s1.undo()
s1.CircleByCenterPerimeter(center=(-5.0, -5.0), point1=(-4.667, 0.0))
s1.undo()
session.viewports['Viewport: 1'].view.setValues(nearPlane=187.116, 
    farPlane=190.008, width=16.0575, height=8.00869, cameraPosition=(-4.09125, 
    -4.95742, 188.562), cameraTarget=(-4.09125, -4.95742, 0))
s1.CircleByCenterPerimeter(center=(-5.0, -5.0), point1=(-4.667, -5.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=187.837, 
    farPlane=189.287, width=6.98231, height=3.48243, cameraPosition=(-5.25126, 
    -4.98086, 188.562), cameraTarget=(-5.25126, -4.98086, 0))
s1.autoTrimCurve(curve1=g[6], point1=(-5.32629346847534, -4.92246341705322))
s1.autoTrimCurve(curve1=g[2], point1=(-5.00410079956055, -4.83161735534668))
s1.autoTrimCurve(curve1=g[5], point1=(-4.79224872589111, -5.00033140182495))
session.viewports['Viewport: 1'].view.setValues(nearPlane=187.709, 
    farPlane=189.415, width=8.21448, height=4.09697, cameraPosition=(-4.62844, 
    -4.78542, 188.562), cameraTarget=(-4.62844, -4.78542, 0))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShell(sketch=s1)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].parts['Part-1'].setValues(geometryRefinement=FINE)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.48, 
    farPlane=29.0886, width=10.6771, height=5.3252, viewOffsetX=-2.04112, 
    viewOffsetY=-2.81726)
p = mdb.models['Model-1'].parts['Part-1']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
    0.004235, 0.004235, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=28.28, gridSpacing=0.7, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.0986, 
    farPlane=30.4699, width=29.1386, height=14.5329, cameraPosition=(1.45946, 
    -0.207906, 28.2843), cameraTarget=(1.45946, -0.207906, 0))
s.CircleByCenterPerimeter(center=(-5.004235, -5.004235), point1=(-3.668, -5.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.4764, 
    farPlane=29.0921, width=10.7698, height=5.37144, cameraPosition=(-2.56504, 
    -3.94181, 28.2843), cameraTarget=(-2.56504, -3.94181, 0))
s.undo()
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2675, 
    farPlane=28.301, width=0.161225, height=0.0804112, cameraPosition=(
    -4.94616, -3.69933, 28.2843), cameraTarget=(-4.94616, -3.69933, 0))
s.ArcByCenterEnds(center=(-5.004235, -5.004235), point1=(-3.675, -5.004235), 
    point2=(-5.004235, -3.675), direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[6], entity2=g[4], addUndoState=False)
s.CoincidentConstraint(entity1=v[7], entity2=g[2], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.6674, 
    farPlane=28.9012, width=8.22467, height=4.10205, cameraPosition=(-2.20287, 
    -4.56124, 28.2843), cameraTarget=(-2.20287, -4.56124, 0))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.7772, 
    farPlane=29.7914, width=19.9309, height=9.94052, viewOffsetX=0.46204, 
    viewOffsetY=-2.27596)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.0335, 
    farPlane=29.535, width=12.1101, height=6.05521, viewOffsetX=-1.32416, 
    viewOffsetY=-3.31995)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#40 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=32, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#1 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=16, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#80 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=17, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#20 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=17, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=SWEEP)
elemType1 = mesh.ElemType(elemCode=CPE8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#2 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.885, 
    farPlane=28.6836, width=3.84591, height=1.92301, viewOffsetX=-3.22797, 
    viewOffsetY=-4.14398)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f1, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[1], sketchPlaneSide=SIDE1, origin=(
    -4.407542, -4.407542, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=28.28, gridSpacing=0.7, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0267, 
    farPlane=28.5419, width=3.42111, height=1.7106, viewOffsetX=-3.71, 
    viewOffsetY=-4.30355)
mdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f, e1 = p.faces, p.edges
p.setSweepPath(region=f[1], edge=e1[7], sense=FORWARD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges1 = e.getSequenceFromMask(mask=('[#20 ]', ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#80 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
    end2Edges=pickedEdges2, ratio=10.0, number=17, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1493, 
    farPlane=28.4192, width=1.79366, height=0.896854, viewOffsetX=-4.15987, 
    viewOffsetY=-4.6487)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges1 = e.getSequenceFromMask(mask=('[#20 ]', ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#80 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
    end2Edges=pickedEdges2, ratio=12.0, number=17, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
#: 
#: Point 1: -3.670765, -5., 0.  Point 2: -3.825081, -5., 0.
#:    Distance: 154.316E-03  Components: -154.316E-03, 0., 0.
#: 
#: Point 1: -3.826496, -4.942349, 0.  Point 2: -3.825081, -5., 0.
#:    Distance: 57.668E-03  Components: 1.415E-03, -57.651E-03, 0.
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f1, e, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[1], sketchPlaneSide=SIDE1, origin=(
    -4.407542, -4.407542, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=28.28, gridSpacing=0.7, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.72452, 
    farPlane=3.79477, width=0.338341, height=0.168747, cameraPosition=(
    -4.90462, -3.89339, 3.75964), cameraTarget=(-4.90462, -3.89339, 0))
s.ArcByCenterEnds(center=(-0.592458, -0.592458), point1=(0.591433534805298, 
    -0.592458), point2=(-0.592458, 0.591821918670655), 
    direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[8], entity2=g[10], addUndoState=False)
s.CoincidentConstraint(entity1=v[9], entity2=g[8], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.41786, 
    farPlane=4.10143, width=4.55671, height=2.27266, cameraPosition=(-3.90183, 
    -4.632, 3.75964), cameraTarget=(-3.90183, -4.632, 0))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#2 ]', ), )
e1, d1 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0606, 
    farPlane=28.5079, width=2.97088, height=1.48548, viewOffsetX=-3.9721, 
    viewOffsetY=-4.45074)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#1 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=32, constraint=FINER)
elemType1 = mesh.ElemType(elemCode=CPE8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#4 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=SWEEP)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0612, 
    farPlane=28.5074, width=2.14651, height=1.07328, viewOffsetX=-3.99423, 
    viewOffsetY=-4.61382)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f, e, d2 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
    -4.199155, -4.199155, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=3.75, gridSpacing=0.09, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9436, 
    farPlane=28.6249, width=4.52193, height=2.26103, viewOffsetX=-3.86536, 
    viewOffsetY=-4.32042)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9943, 
    farPlane=28.5743, width=3.85061, height=1.92536, viewOffsetX=-3.99618, 
    viewOffsetY=-4.20466)
p = mdb.models['Model-1'].parts['Part-1']
p.deleteMesh()
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#4 ]', ), )
p.PartitionEdgeByParam(edges=pickedEdges, parameter=0.915391519370898)
p = mdb.models['Model-1'].parts['Part-1']
del p.features['Partition edge-1']
p = mdb.models['Model-1'].parts['Part-1']
f1 = p.faces
e1 = p.edges
p.PartitionFaceByCurvedPathEdgeParams(face=f1[0], edge1=e1[0], edge2=e1[2], 
    parameter1=0.5, parameter2=1.0)
#* Feature creation failed.
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['Part-1']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
    -4.199155, -4.199155, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=3.75, 
    gridSpacing=0.09, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7138, 
    farPlane=3.80548, width=0.441544, height=0.22022, cameraPosition=(-4.90851, 
    -3.8385, 3.75964), cameraTarget=(-4.90851, -3.8385, 0))
s.ArcByCenterEnds(center=(-0.800845, -0.800845), point1=(0.452872422485352, 
    -0.800845), point2=(-0.800845, 0.45), direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[10], entity2=g[3], addUndoState=False)
s.CoincidentConstraint(entity1=v[11], entity2=g[5], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.38809, 
    farPlane=4.1312, width=4.95355, height=2.47058, cameraPosition=(-3.93812, 
    -4.3238, 3.75964), cameraTarget=(-3.93812, -4.3238, 0))
s.radialPattern(geomList=(g[5], ), vertexList=(), number=32, totalAngle=-90.0, 
    centerPoint=(-0.800845, -0.800845))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.61951, 
    farPlane=3.89978, width=1.34979, height=0.673211, cameraPosition=(-3.79808, 
    -4.81138, 3.75964), cameraTarget=(-3.79808, -4.81138, 0))
s.Line(point1=(0.452131131965858, -0.757738293245818), point2=(0.495, -0.7875))
s.CoincidentConstraint(entity1=v[74], entity2=g[15], addUndoState=False)
s.undo()
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.66988, 
    farPlane=3.84941, width=0.996917, height=0.497213, cameraPosition=(
    -3.78479, -4.84586, 3.75964), cameraTarget=(-3.78479, -4.84586, 0))
s.Line(point1=(0.451263285164066, -0.737345254592504), point2=(0.52839, 
    -0.800845))
s.CoincidentConstraint(entity1=v[74], entity2=g[15], addUndoState=False)
s.CoincidentConstraint(entity1=v[73], entity2=g[3], addUndoState=False)
s.Line(point1=(0.451263285164066, -0.737345254592504), point2=(
    0.521570125078826, -0.666368525523036))
s.undo()
s.Line(point1=(0.446440003832633, -0.674008512120749), point2=(
    0.513061693088277, -0.599561898168704))
s.CoincidentConstraint(entity1=v[75], entity2=g[15], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.40963, 
    farPlane=4.10966, width=3.37149, height=1.68153, cameraPosition=(-3.88552, 
    -4.53899, 3.75964), cameraTarget=(-3.88552, -4.53899, 0))
s.undo()
s.undo()
s.undo()
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.44986, 
    farPlane=4.06943, width=3.58215, height=1.7866, cameraPosition=(-3.85201, 
    -4.54401, 3.75964), cameraTarget=(-3.85201, -4.54401, 0))
s.radialPattern(geomList=(g[5], ), vertexList=(), number=33, totalAngle=-90.0, 
    centerPoint=(0.0, 0.0))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.24492, 
    farPlane=4.27437, width=6.86228, height=3.42256, cameraPosition=(-3.06947, 
    -4.81468, 3.75964), cameraTarget=(-3.06947, -4.81468, 0))
s.undo()
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.49095, 
    farPlane=4.02833, width=3.58215, height=1.7866, cameraPosition=(-3.86775, 
    -4.30274, 3.75964), cameraTarget=(-3.86775, -4.30274, 0))
s.radialPattern(geomList=(g[5], ), vertexList=(), number=33, totalAngle=-90.0, 
    centerPoint=(-0.800845, -0.800845))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.53126, 
    farPlane=3.98803, width=2.19989, height=1.0972, cameraPosition=(-3.87249, 
    -4.57955, 3.75964), cameraTarget=(-3.87249, -4.57955, 0))
s.Line(point1=(0.45136226494363, -0.739328001814879), point2=(0.52839, 
    -0.800845))
s.CoincidentConstraint(entity1=v[76], entity2=g[15], addUndoState=False)
s.CoincidentConstraint(entity1=v[75], entity2=g[3], addUndoState=False)
s.Line(point1=(0.439302824456101, -0.616886447765771), point2=(
    0.502849122196788, -0.541524115814892))
s.CoincidentConstraint(entity1=v[77], entity2=g[15], addUndoState=False)
s.autoTrimCurve(curve1=g[44], point1=(0.469895270614624, -0.607327791900635))
s.autoTrimCurve(curve1=g[46], point1=(0.504659560470581, -0.738173338623047))
s.autoTrimCurve(curve1=g[15], point1=(0.435130742340088, -0.590972277374267))
s.autoTrimCurve(curve1=g[52], point1=(0.450427201538086, -0.779063078613281))
s.Line(point1=(0.415300082485299, -0.496216515136653), point2=(
    0.502849122196788, -0.541524115814892))
s.CoincidentConstraint(entity1=v[83], entity2=g[42], addUndoState=False)
s.Line(point1=(0.379585198557637, -0.378480321343638), point2=(
    0.427208010397641, -0.29216878778019))
s.CoincidentConstraint(entity1=v[84], entity2=g[40], addUndoState=False)
s.Line(point1=(0.335074739864205, -0.2635949748046), point2=(0.427208010397641, 
    -0.29216878778019))
s.CoincidentConstraint(entity1=v[85], entity2=g[38], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[57], entity2=v[56], midpoint=v[85], 
    addUndoState=False)
s.Line(point1=(0.274504302121485, -0.156305432657428), point2=(
    0.304373510108974, -0.062361601312189))
s.CoincidentConstraint(entity1=v[86], entity2=g[36], addUndoState=False)
s.Line(point1=(0.304373510108974, -0.062361601312189), point2=(
    0.206150276088736, -0.0540064033954298))
s.CoincidentConstraint(entity1=v[87], entity2=g[34], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.51699, 
    farPlane=4.0023, width=2.33734, height=1.16575, cameraPosition=(-4.11187, 
    -4.3662, 3.75964), cameraTarget=(-4.11187, -4.3662, 0))
s.Line(point1=(0.12809833506764, 0.0411001619177606), point2=(
    0.139066082290501, 0.139066082290501))
s.CoincidentConstraint(entity1=v[88], entity2=g[32], addUndoState=False)
s.Line(point1=(0.139066082290501, 0.139066082290501), point2=(
    0.0411001619177609, 0.12809833506764))
s.CoincidentConstraint(entity1=v[89], entity2=g[30], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(width=2.74982, height=1.37147, 
    cameraPosition=(-4.0871, -4.41277, 3.75964), cameraTarget=(-4.0871, 
    -4.41277, 0))
s.Line(point1=(-0.0540064033954295, 0.206150276088736), point2=(
    -0.0623616013121887, 0.304373510108974))
s.CoincidentConstraint(entity1=v[90], entity2=g[28], addUndoState=False)
s.Line(point1=(-0.0623616013121887, 0.304373510108974), point2=(
    -0.156305432657428, 0.274504302121485))
s.CoincidentConstraint(entity1=v[91], entity2=g[26], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.61062, 
    farPlane=3.90867, width=1.98674, height=0.990889, cameraPosition=(-4.26533, 
    -4.19809, 3.75964), cameraTarget=(-4.26533, -4.19809, 0))
s.Line(point1=(-0.264811730294102, 0.332502126529079), point2=(
    -0.292168787780189, 0.42720801039764))
s.CoincidentConstraint(entity1=v[92], entity2=g[24], addUndoState=False)
s.Line(point1=(-0.292168787780189, 0.42720801039764), point2=(
    -0.378480321343637, 0.379585198557637))
s.CoincidentConstraint(entity1=v[93], entity2=g[22], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.55339, 
    farPlane=3.9659, cameraPosition=(-4.31266, -4.14259, 3.75964), 
    cameraTarget=(-4.31266, -4.14259, 0))
s.Line(point1=(-0.496216515136653, 0.4153000824853), point2=(
    -0.541524115814892, 0.502849122196788))
s.CoincidentConstraint(entity1=v[94], entity2=g[20], addUndoState=False)
s.Line(point1=(-0.541524115814892, 0.502849122196788), point2=(
    -0.61646887559083, 0.442117867399294))
s.CoincidentConstraint(entity1=v[95], entity2=g[18], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[17], entity2=v[16], midpoint=v[95], 
    addUndoState=False)
s.Line(point1=(-0.739328001814879, 0.45136226494363), point2=(-0.800845, 
    0.52839))
s.CoincidentConstraint(entity1=v[96], entity2=g[16], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.55339, 
    farPlane=3.9659, cameraPosition=(-3.72946, -4.60062, 3.75964), 
    cameraTarget=(-3.72946, -4.60062, 0))
s.autoTrimCurve(curve1=g[53], point1=(0.423231747894287, -0.536251875610351))
s.autoTrimCurve(curve1=g[69], point1=(0.367974665908814, -0.351613852233887))
s.autoTrimCurve(curve1=g[41], point1=(0.42699923828125, -0.429161879272461))
s.undo()
s.autoTrimCurve(curve1=g[71], point1=(0.341601994781494, -0.285144659729004))
s.autoTrimCurve(curve1=g[42], point1=(0.477232840805054, -0.488245817871094))
s.autoTrimCurve(curve1=g[40], point1=(0.41192903831482, -0.3639234034729))
session.viewports['Viewport: 1'].view.setValues(width=2.33734, height=1.16575, 
    cameraPosition=(-3.69555, -4.58227, 3.75964), cameraTarget=(-3.69555, 
    -4.58227, 0))
s.autoTrimCurve(curve1=g[38], point1=(0.361769345550537, -0.247710081787109))
s.autoTrimCurve(curve1=g[36], point1=(0.304148343353272, -0.140547606201172))
s.autoTrimCurve(curve1=g[34], point1=(0.250959780960083, -0.0319374530029295))
s.autoTrimCurve(curve1=g[72], point1=(0.252437260894776, -0.124618384094238))
s.autoTrimCurve(curve1=g[79], point1=(0.230275300292969, -0.0898626773071287))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58433, 
    farPlane=3.93496, width=1.68873, height=0.842255, cameraPosition=(-4.063, 
    -4.28162, 3.75964), cameraTarget=(-4.063, -4.28162, 0))
s.autoTrimCurve(curve1=g[32], point1=(0.158571389465332, 0.0718599827575686))
s.autoTrimCurve(curve1=g[30], point1=(0.062499192504883, 0.149284508972168))
s.autoTrimCurve(curve1=g[80], point1=(0.118007329254151, 0.0645357640075686))
s.autoTrimCurve(curve1=g[84], point1=(0.0678364308166506, 0.0969707043457033))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58433, 
    farPlane=3.93496, width=2.33734, height=1.16575, cameraPosition=(-4.11183, 
    -4.32006, 3.75964), cameraTarget=(-4.11183, -4.32006, 0))
s.autoTrimCurve(curve1=g[28], point1=(-0.0353086917114256, 0.237509396820069))
s.autoTrimCurve(curve1=g[26], point1=(-0.135776373596191, 0.30702009513855))
s.autoTrimCurve(curve1=g[85], point1=(-0.0796326129150389, 0.221579936294556))
s.autoTrimCurve(curve1=g[89], point1=(-0.128389212341308, 0.25778332069397))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58433, 
    farPlane=3.93496, width=1.68873, height=0.842255, cameraPosition=(-4.31336, 
    -4.09587, 3.75964), cameraTarget=(-4.31336, -4.09587, 0))
s.autoTrimCurve(curve1=g[24], point1=(-0.241233679504394, 0.377933409957886))
s.autoTrimCurve(curve1=g[22], point1=(-0.359722468109131, 0.420831111221314))
s.autoTrimCurve(curve1=g[90], point1=(-0.339440676422119, 0.361193087844849))
s.autoTrimCurve(curve1=g[93], point1=(-0.292472216339111, 0.351776507644654))
s.autoTrimCurve(curve1=g[20], point1=(-0.48568329498291, 0.450126794128418))
s.autoTrimCurve(curve1=g[18], point1=(-0.609509798736572, 0.472098735122681))
s.autoTrimCurve(curve1=g[94], point1=(-0.588160368652344, 0.43024745300293))
s.autoTrimCurve(curve1=g[98], point1=(-0.539056631774902, 0.429201272277832))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58433, 
    farPlane=3.93496, height=0.842255, cameraPosition=(-4.34976, -4.04971, 
    3.75964), cameraTarget=(-4.34976, -4.04971, 0))
s.autoTrimCurve(curve1=g[16], point1=(-0.737712713928222, 0.497332719116211))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.46499, 
    farPlane=4.0543, width=3.92831, height=1.95925, cameraPosition=(-4.18271, 
    -4.06616, 3.75964), cameraTarget=(-4.18271, -4.06616, 0))
mdb.saveAs(
    pathName='/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen')
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-1']
s1 = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s1)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=20.1147, 
    farPlane=20.2371, width=0.589553, height=0.29404, cameraPosition=(-4.94182, 
    -3.67847, 20.1759), cameraTarget=(-4.94182, -3.67847, 0))
s2.autoTrimCurve(curve1=g[99], point1=(-0.755336138458252, 0.451099065093994))
s2.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s2)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#2 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #20000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#10 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#8 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1744, 
    farPlane=28.3941, width=1.05601, height=0.528019, viewOffsetX=-4.58076, 
    viewOffsetY=-3.84316)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#8 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#8 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1323, 
    farPlane=28.4363, width=2.01995, height=1.01, viewOffsetX=-4.34773, 
    viewOffsetY=-3.98972)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#7fe0 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0374, 
    farPlane=28.5311, width=3.27806, height=1.63908, viewOffsetX=-4.15539, 
    viewOffsetY=-4.09816)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#8 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#ffff8000 #1 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(width=3.85061, height=1.92536, 
    viewOffsetX=-4.03772, viewOffsetY=-4.06691)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #fffe ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9943, 
    farPlane=28.5743, width=3.85061, height=1.92536, viewOffsetX=-4.02873, 
    viewOffsetY=-4.05011)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1594, 
    farPlane=28.4092, width=1.20317, height=0.600081, cameraPosition=(-3.67295, 
    -4.47263, 28.2843), cameraTarget=(-3.67295, -4.47263, 0))
s1.autoTrimCurve(curve1=g[57], point1=(0.352802422790528, -0.268630835266113))
s1.autoTrimCurve(curve1=g[39], point1=(0.359647420196533, -0.314848753662109))
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s1)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
#: Warning: Mesh deleted in 4 regions due to geometry association failure.
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1463, 
    farPlane=28.4222, width=1.3287, height=0.662691, cameraPosition=(-4.39666, 
    -3.67148, 28.2843), cameraTarget=(-4.39666, -3.67148, 0))
s2.autoTrimCurve(curve1=g[67], point1=(-0.609894129486084, 0.453994420318604))
s2.autoTrimCurve(curve1=g[66], point1=(-0.513306948394775, 0.433414128570557))
s2.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s2)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
#: Warning: Mesh deleted in 42 regions due to geometry association failure.
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.02, 
    farPlane=28.5485, width=2.54538, height=1.26951, cameraPosition=(-3.26604, 
    -4.47185, 28.2843), cameraTarget=(-3.26604, -4.47185, 0))
s1.Line(point1=(0.332502126529078, -0.264811730294103), point2=(
    0.427208010397641, -0.29216878778019))
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2322, 
    farPlane=28.3363, width=0.50112, height=0.249934, cameraPosition=(-3.76224, 
    -4.52948, 28.2843), cameraTarget=(-3.76224, -4.52948, 0))
s1.Line(point1=(0.427208010397641, -0.29216878778019), point2=(
    0.34378710105896, -0.331721159667969))
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9185, 
    farPlane=28.65, width=4.87615, height=2.43198, cameraPosition=(-2.99357, 
    -4.27913, 28.2843), cameraTarget=(-2.99357, -4.27913, 0))
s1.Line(point1=(0.330166247634888, -0.35469660446167), point2=(0.1125, -0.45))
s1.autoTrimCurve(curve1=g[68], point1=(-0.773228022308349, 0.444180634765625))
s1.autoTrimCurve(curve1=g[101], point1=(-0.745488020629883, 0.41396966293335))
s1.autoTrimCurve(curve1=g[17], point1=(-0.671513411254883, 0.465328362731934))
s1.autoTrimCurve(curve1=g[141], point1=(-0.674595686645508, 0.416990664749146))
s1.autoTrimCurve(curve1=g[110], point1=(-0.702336165161133, 0.444180634765625))
s1.autoTrimCurve(curve1=g[97], point1=(-0.612950178833008, 0.407927659301758))
s1.autoTrimCurve(curve1=g[19], point1=(-0.557469221801758, 0.407927659301758))
s1.autoTrimCurve(curve1=g[96], point1=(-0.508152815551758, 0.386779931335449))
s1.autoTrimCurve(curve1=g[100], point1=(-0.480412813873291, 0.410948661117554))
s1.autoTrimCurve(curve1=g[21], point1=(-0.446507784576416, 0.359590199737549))
s1.autoTrimCurve(curve1=g[143], point1=(-0.443425509185791, 0.435117390899658))
s1.autoTrimCurve(curve1=g[92], point1=(-0.378698202819824, 0.347505715637207))
s1.autoTrimCurve(curve1=g[65], point1=(-0.363286825866699, 0.377716687469483))
s1.autoTrimCurve(curve1=g[142], point1=(-0.400274130554199, 0.389801171569824))
s1.autoTrimCurve(curve1=g[64], point1=(-0.270818564147949, 0.368653205184937))
s1.autoTrimCurve(curve1=g[91], point1=(-0.289312216491699, 0.33240022972107))
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
del p.features['Partition face-3']
#: Warning: Failed to attach mesh to part geometry.
p = mdb.models['Model-1'].parts['Part-1']
f1, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[0], sketchPlaneSide=SIDE1, origin=(
    -4.199155, -4.199155, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=3.75, 
    gridSpacing=0.09, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.73427, 
    farPlane=3.78502, width=0.244451, height=0.12192, cameraPosition=(-4.9354, 
    -3.77321, 3.75964), cameraTarget=(-4.9354, -3.77321, 0))
s.ArcByCenterEnds(center=(-0.800845, -0.800845), point1=(0.455016520767212, 
    -0.800845), point2=(-0.800845, 0.454952147750855), 
    direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[10], entity2=g[3], addUndoState=False)
s.CoincidentConstraint(entity1=v[11], entity2=g[5], addUndoState=False)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-1']
s1 = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s1)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.68048, 
    farPlane=3.83881, width=0.762535, height=0.380314, cameraPosition=(
    -4.88088, -3.82994, 3.75964), cameraTarget=(-4.88088, -3.82994, 0))
s2.Line(point1=(-0.800845, 0.52839), point2=(-0.800845, 0.383046534805298))
s2.VerticalConstraint(entity=g[27], addUndoState=False)
s2.PerpendicularConstraint(entity1=g[4], entity2=g[27], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.41786, 
    farPlane=4.10143, width=4.55671, height=2.27266, cameraPosition=(-3.9607, 
    -4.24359, 3.75964), cameraTarget=(-3.9607, -4.24359, 0))
s2.radialPattern(geomList=(g[5], ), vertexList=(), number=33, totalAngle=-90.0, 
    centerPoint=(-0.800845, -0.800845))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7238, 
    farPlane=3.79549, width=0.345246, height=0.172191, cameraPosition=(
    -3.71636, -4.93552, 3.75964), cameraTarget=(-3.71636, -4.93552, 0))
s2.Line(point1=(0.455016520767212, -0.800845), point2=(0.526788878228883, 
    -0.735622529915394))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.67886, 
    farPlane=3.84042, width=1.07695, height=0.537129, cameraPosition=(-3.77534, 
    -4.85487, 3.75964), cameraTarget=(-3.77534, -4.85487, 0))
s2.autoTrimCurve(curve1=g[60], point1=(0.483037141113281, -0.766809317321777))
s2.Line(point1=(0.453503780565209, -0.739222795898658), point2=(0.52839, 
    -0.800845))
s2.CoincidentConstraint(entity1=v[79], entity2=g[15], addUndoState=False)
s2.CoincidentConstraint(entity1=v[78], entity2=g[3], addUndoState=False)
s2.Line(point1=(0.441423716111573, -0.616571843207595), point2=(
    0.502849122196788, -0.541524115814892))
s2.CoincidentConstraint(entity1=v[80], entity2=g[15], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.72224, 
    farPlane=3.79705, width=0.498681, height=0.248717, cameraPosition=(
    -3.75198, -4.99067, 3.75964), cameraTarget=(-3.75198, -4.99067, 0))
s2.autoTrimCurve(curve1=g[61], point1=(0.479956773071289, -0.757061335296631))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7586, 
    farPlane=3.76069, width=0.0100897, height=0.00503223, cameraPosition=(
    -3.74635, -4.93927, 3.75964), cameraTarget=(-3.74635, -4.93927, 0))
s2.Line(point1=(0.52839, -0.800845), point2=(0.453379239791078, 
    -0.739228914194559))
s2.CoincidentConstraint(entity1=v[78], entity2=g[3], addUndoState=False)
s2.CoincidentConstraint(entity1=v[81], entity2=g[58], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.62575, 
    farPlane=3.89354, width=1.2897, height=0.643236, cameraPosition=(-3.87132, 
    -4.95558, 3.75964), cameraTarget=(-3.87132, -4.95558, 0))
s2.autoTrimCurve(curve1=g[62], point1=(0.481917527465821, -0.571440550537109))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.72199, 
    farPlane=3.7973, width=0.502043, height=0.250394, cameraPosition=(-3.78166, 
    -4.82943, 3.75964), cameraTarget=(-3.78166, -4.82943, 0))
s2.Line(point1=(0.441143987681443, -0.616613337000521), point2=(
    0.502849122196788, -0.541524115814892))
s2.CoincidentConstraint(entity1=v[82], entity2=g[56], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7584, 
    farPlane=3.76089, width=0.0119502, height=0.00596018, cameraPosition=(
    -3.78146, -4.69691, 3.75964), cameraTarget=(-3.78146, -4.69691, 0))
s2.Line(point1=(0.502849122196788, -0.541524115814892), point2=(
    0.417106287763241, -0.4957640842671))
s2.CoincidentConstraint(entity1=v[83], entity2=g[54], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.68032, 
    farPlane=3.83897, width=1.05751, height=0.527433, cameraPosition=(-3.9452, 
    -4.65775, 3.75964), cameraTarget=(-3.9452, -4.65775, 0))
s2.Line(point1=(0.381483568820775, -0.377801073602314), point2=(
    0.427208010397641, -0.29216878778019))
s2.CoincidentConstraint(entity1=v[84], entity2=g[52], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7595, 
    farPlane=3.75978, width=0.0013504, height=0.00067351, cameraPosition=(
    -3.86486, -4.4631, 3.75964), cameraTarget=(-3.86486, -4.4631, 0))
s2.Line(point1=(0.427208010397641, -0.29216878778019), point2=(
    0.33432224801677, -0.263950876942582))
s2.CoincidentConstraint(entity1=v[85], entity2=g[50], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.67419, 
    farPlane=3.8451, width=1.13934, height=0.568245, cameraPosition=(-4.14537, 
    -4.49695, 3.75964), cameraTarget=(-4.14537, -4.49695, 0))
s2.Line(point1=(0.276103125991921, -0.155347134508468), point2=(
    0.304373510108974, -0.062361601312189))
s2.CoincidentConstraint(entity1=v[86], entity2=g[48], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75961, 
    farPlane=3.75967, width=0.00028643, height=0.000142857, cameraPosition=(
    -3.99155, -4.25206, 3.75964), cameraTarget=(-3.99155, -4.25206, 0))
s2.Line(point1=(0.304373510108974, -0.062361601312189), point2=(
    0.207647384390016, -0.0528960722059577))
s2.CoincidentConstraint(entity1=v[87], entity2=g[46], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.56618, 
    farPlane=3.95311, width=2.57923, height=1.28639, cameraPosition=(-4.59866, 
    -4.38306, 3.75964), cameraTarget=(-4.59866, -4.38306, 0))
s2.Line(point1=(0.12959584992461, 0.0424574302689395), point2=(
    0.139066082290501, 0.139066082290501))
s2.CoincidentConstraint(entity1=v[88], entity2=g[44], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75874, 
    farPlane=3.76055, width=0.00873281, height=0.00435549, cameraPosition=(
    -4.15746, -4.06917, 3.75964), cameraTarget=(-4.15746, -4.06917, 0))
s2.Line(point1=(0.139066082290501, 0.139066082290501), point2=(
    0.0424501015302398, 0.129587763907518))
s2.CoincidentConstraint(entity1=v[89], entity2=g[42], addUndoState=False)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75504, 
    farPlane=3.76425, width=0.0443572, height=0.0221232, cameraPosition=(
    -4.25235, -3.99348, 3.75964), cameraTarget=(-4.25235, -3.99348, 0))
s2.Line(point1=(-0.052729165540162, 0.207872431976961), point2=(
    -0.0496867625427244, 0.20708694770813))
s2.CoincidentConstraint(entity1=v[90], entity2=g[15], addUndoState=False)
s2.autoTrimCurve(curve1=g[72], point1=(-0.0504716365051268, 0.207526829986572))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.70634, 
    farPlane=3.81295, width=0.710616, height=0.35442, cameraPosition=(-4.3069, 
    -3.9796, 3.75964), cameraTarget=(-4.3069, -3.9796, 0))
s2.Line(point1=(-0.0528959733693601, 0.20764751765574), point2=(
    -0.0623616013121887, 0.304373510108974))
s2.CoincidentConstraint(entity1=v[92], entity2=g[40], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75943, 
    farPlane=3.75986, width=0.00204511, height=0.00102, cameraPosition=(
    -4.35455, -3.92275, 3.75964), cameraTarget=(-4.35455, -3.92275, 0))
s2.Line(point1=(-0.0623616013121887, 0.304373510108974), point2=(
    -0.155346988894097, 0.276103368934822))
s2.CoincidentConstraint(entity1=v[93], entity2=g[38], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.687, 
    farPlane=3.83228, width=0.968437, height=0.483008, cameraPosition=(
    -4.36469, -3.85815, 3.75964), cameraTarget=(-4.36469, -3.85815, 0))
s2.Line(point1=(-0.26395099652482, 0.33432199518137), point2=(
    -0.292168787780189, 0.42720801039764))
s2.CoincidentConstraint(entity1=v[94], entity2=g[36], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75924, 
    farPlane=3.76005, width=0.00533922, height=0.00266294, cameraPosition=(
    -4.57691, -3.818, 3.75964), cameraTarget=(-4.57691, -3.818, 0))
s2.Line(point1=(-0.292168787780189, 0.42720801039764), point2=(
    -0.377802567413971, 0.381479393896875))
s2.CoincidentConstraint(entity1=v[95], entity2=g[34], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.69917, 
    farPlane=3.82012, width=0.8063, height=0.402142, cameraPosition=(-4.55726, 
    -3.78671, 3.75964), cameraTarget=(-4.55726, -3.78671, 0))
s2.Line(point1=(-0.495763436986236, 0.417108871853299), point2=(
    -0.541524115814892, 0.502849122196788))
s2.CoincidentConstraint(entity1=v[96], entity2=g[32], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.7595, 
    farPlane=3.75979, width=0.00142507, height=0.000710753, cameraPosition=(
    -4.81534, -3.75801, 3.75964), cameraTarget=(-4.81534, -3.75801, 0))
s2.Line(point1=(-0.541524115814892, 0.502849122196788), point2=(
    -0.616612849530385, 0.441147273938161))
s2.CoincidentConstraint(entity1=v[97], entity2=g[30], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.69868, 
    farPlane=3.82061, width=0.812789, height=0.405378, cameraPosition=(
    -4.77591, -3.81503, 3.75964), cameraTarget=(-4.77591, -3.81503, 0))
s2.Line(point1=(-0.739229258659055, 0.453372228055186), point2=(-0.800845, 
    0.52839))
s2.CoincidentConstraint(entity1=v[98], entity2=g[28], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.598, 
    farPlane=3.92129, width=1.55705, height=0.776578, cameraPosition=(-3.83529, 
    -4.72178, 3.75964), cameraTarget=(-3.83529, -4.72178, 0))
s2.autoTrimCurve(curve1=g[15], point1=(0.457363751678467, -0.774895521850586))
s2.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s2)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9943, 
    farPlane=28.5743, width=3.85061, height=1.92536, viewOffsetX=-3.87945, 
    viewOffsetY=-4.28746)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=20.0136, 
    farPlane=20.3382, width=1.56318, height=0.779637, cameraPosition=(-3.86219, 
    -4.7272, 20.1759), cameraTarget=(-3.86219, -4.7272, 0))
s1.autoTrimCurve(curve1=g[56], point1=(0.475295213012696, -0.615697714538574))
s1.autoTrimCurve(curve1=g[54], point1=(0.454545167236328, -0.492699476928711))
s1.autoTrimCurve(curve1=g[80], point1=(0.427866366653443, -0.535312983245849))
s1.autoTrimCurve(curve1=g[91], point1=(0.436759141235352, -0.586643072814941))
session.viewports['Viewport: 1'].view.setValues(nearPlane=19.985, 
    farPlane=20.3668, width=2.54538, height=1.26951, cameraPosition=(-3.94329, 
    -4.72359, 20.1759), cameraTarget=(-3.94329, -4.72359, 0))
s1.autoTrimCurve(curve1=g[52], point1=(0.408721116333008, -0.365944239349365))
s1.autoTrimCurve(curve1=g[50], point1=(0.38136878326416, -0.244513365478515))
s1.autoTrimCurve(curve1=g[92], point1=(0.354016450195313, -0.296555372924804))
s1.autoTrimCurve(curve1=g[96], point1=(0.370105889587403, -0.343866202087402))
session.viewports['Viewport: 1'].view.setValues(nearPlane=19.985, 
    farPlane=20.3668, width=1.83904, height=0.91722, cameraPosition=(-3.95858, 
    -4.46794, 20.1759), cameraTarget=(-3.95858, -4.46794, 0))
s1.autoTrimCurve(curve1=g[48], point1=(0.298698571472168, -0.141745897979736))
s1.autoTrimCurve(curve1=g[46], point1=(0.232437279968262, -0.0300844638061522))
s1.autoTrimCurve(curve1=g[97], point1=(0.226625111846924, -0.079079004974365))
s1.autoTrimCurve(curve1=g[101], point1=(0.252199557571411, -0.126933905334472))
session.viewports['Viewport: 1'].view.setValues(nearPlane=20.038, 
    farPlane=20.3139, width=1.3287, height=0.662691, cameraPosition=(-4.1133, 
    -4.24585, 20.1759), cameraTarget=(-4.1133, -4.24585, 0))
s1.autoTrimCurve(curve1=g[42], point1=(0.0657025845336916, 0.147995618133545))
s1.autoTrimCurve(curve1=g[44], point1=(0.17068877532959, 0.0739060910034182))
s1.autoTrimCurve(curve1=g[102], point1=(0.107697156219483, 0.0673204930114748))
s1.autoTrimCurve(curve1=g[107], point1=(0.0766212017822268, 
    0.0969563992309572))
session.viewports['Viewport: 1'].view.setValues(nearPlane=20.0136, 
    farPlane=20.3382, width=1.56318, height=0.779637, cameraPosition=(-4.16771, 
    -4.19059, 20.1759), cameraTarget=(-4.16771, -4.19059, 0))
s1.autoTrimCurve(curve1=g[38], point1=(-0.142456385345459, 0.299595263748169))
s1.autoTrimCurve(curve1=g[40], point1=(-0.0317886798095701, 0.237611678390503))
s1.autoTrimCurve(curve1=g[108], point1=(-0.0703247515869139, 
    0.223084357528687))
s1.autoTrimCurve(curve1=g[112], point1=(-0.127634855957031, 0.25407615020752))
session.viewports['Viewport: 1'].view.setValues(nearPlane=19.9513, 
    farPlane=20.4005, width=2.99456, height=1.49354, cameraPosition=(-4.01208, 
    -4.33825, 20.1759), cameraTarget=(-4.01208, -4.33825, 0))
s1.autoTrimCurve(curve1=g[36], point1=(-0.246398779602051, 0.379469302444458))
s1.autoTrimCurve(curve1=g[34], point1=(-0.363758417816162, 0.411009934692383))
s1.autoTrimCurve(curve1=g[113], point1=(-0.344829413146972, 0.370192673950196))
s1.autoTrimCurve(curve1=g[116], point1=(-0.295614096374512, 0.357205298690796))
s1.autoTrimCurve(curve1=g[32], point1=(-0.490582796783447, 0.446261075286865))
s1.autoTrimCurve(curve1=g[117], point1=(-0.52654776260376, 0.429563191680908))
s1.autoTrimCurve(curve1=g[121], point1=(-0.575763079376221, 0.435129073410034))
s1.autoTrimCurve(curve1=g[30], point1=(-0.607942434997558, 0.479657080917359))
s1.autoTrimCurve(curve1=g[28], point1=(-0.736658903808594, 0.487078336029053))
s1.autoTrimCurve(curve1=g[122], point1=(-0.78208861038208, 0.461103823928833))
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s1)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.8024, 
    farPlane=28.7661, width=4.64369, height=2.32191, viewOffsetX=-3.38664, 
    viewOffsetY=-4.10253)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#bdbdbdbd #fbbdbdbd ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.7998, 
    farPlane=28.7688, width=6.42664, height=3.21341, viewOffsetX=-3.74049, 
    viewOffsetY=-4.10662)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1463, 
    farPlane=28.4222, width=1.3287, height=0.662691, cameraPosition=(-3.49778, 
    -4.74274, 28.2843), cameraTarget=(-3.49778, -4.74274, 0))
s2.autoTrimCurve(curve1=g[58], point1=(0.497278359680176, -0.739925238342285))
s2.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s2)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
#: Warning: Mesh deleted in 2 regions due to geometry association failure.
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9264, 
    farPlane=28.6422, width=4.12723, height=2.06368, viewOffsetX=-4.61253, 
    viewOffsetY=-3.95316)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#dededede #7ddedede ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#dededede #7ddedede ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#dededede #7ddedede ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#dededede #7ddedede ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9185, 
    farPlane=28.65, width=4.85421, height=2.42717, viewOffsetX=-3.79812, 
    viewOffsetY=-4.0791)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0:2 #2 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0192, 
    farPlane=28.5493, width=3.51981, height=1.75996, viewOffsetX=-3.74368, 
    viewOffsetY=-4.48719)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1454, 
    farPlane=28.4232, width=1.84564, height=0.922846, viewOffsetX=-4.34127, 
    viewOffsetY=-4.1158)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1846, 
    farPlane=28.3839, width=0.959989, height=0.478795, cameraPosition=(
    -4.66654, -3.65152, 28.2843), cameraTarget=(-4.66654, -3.65152, 0))
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0587, 
    farPlane=28.5098, width=2.99606, height=1.49807, viewOffsetX=-4.0541, 
    viewOffsetY=-4.1902)
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1800000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:4 #40000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1800000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #d000000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:4 #500000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #d000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1839, 
    farPlane=28.3846, width=0.964754, height=0.482391, viewOffsetX=-4.56859, 
    viewOffsetY=-3.86339)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #5d000000 #2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:4 #6c500000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2116, 
    farPlane=28.3569, width=0.965702, height=0.482865, viewOffsetX=-4.57855, 
    viewOffsetY=-3.86837)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40000000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:5 #2 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2131, 
    farPlane=28.3554, width=0.683802, height=0.341911, viewOffsetX=-4.65134, 
    viewOffsetY=-3.83293)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4000000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:4 #10000000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
p.generateMesh(regions=pickedRegions)
elemType1 = mesh.ElemType(elemCode=CPE8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
p = mdb.models['Model-1'].parts['Part-1']
p.deleteMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2131, 
    farPlane=28.3555, width=0.946437, height=0.473232, viewOffsetX=-4.5736, 
    viewOffsetY=-3.82469)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=FREE)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000000 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #800000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #100000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2131, 
    farPlane=28.3555, width=0.946437, height=0.473232, viewOffsetX=-4.52117, 
    viewOffsetY=-3.82742)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #140000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #400000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.8476, 
    farPlane=28.7209, width=5.79352, height=2.89684, viewOffsetX=-3.97284, 
    viewOffsetY=-3.8145)
p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()
#*p.generateMesh()
#* Command Interrupted
p = mdb.models['Model-1'].parts['Part-1']
p.generateMesh()
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9676, 
    farPlane=28.6009, width=4.20386, height=2.10199, viewOffsetX=-3.70695, 
    viewOffsetY=-4.06521)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#dedede00 #7ddedede ]', ), )
p.generateMesh(regions=pickedRegions)
elemType1 = mesh.ElemType(elemCode=CPE8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#0 #7dde0000 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1432, 
    farPlane=28.4253, width=1.35612, height=0.678082, viewOffsetX=-4.53384, 
    viewOffsetY=-3.85144)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40000000 ]', ), )
p.generateMesh(regions=pickedRegions)
mdb.meshEditOptions.setValues(enableUndo=True, maxUndoCacheElements=0.5)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #7ddedeca ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=(
    '[#0:2 #41000000 #8092140 #a4110124 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #20000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.7182, 
    farPlane=28.8504, width=7.50534, height=3.75278, viewOffsetX=-2.68214, 
    viewOffsetY=-3.79683)
p = mdb.models['Model-1'].parts['Part-1']
del p.features['Partition face-3']
#: Warning: Failed to attach mesh to part geometry.
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.8019, 
    farPlane=28.7667, width=6.39879, height=3.19949, viewOffsetX=-3.08423, 
    viewOffsetY=-3.74822)
p = mdb.models['Model-1'].parts['Part-1']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
    -4.199155, -4.199155, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=3.75, 
    gridSpacing=0.09, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.72653, 
    farPlane=3.79276, width=0.319013, height=0.159108, cameraPosition=(
    -4.95158, -3.75957, 3.75964), cameraTarget=(-4.95158, -3.75957, 0))
s.ArcByCenterEnds(center=(-0.800845, -0.800845), point1=(0.455276635437012, 
    -0.800845), point2=(-0.800845, 0.455718267402649), 
    direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[10], entity2=g[3], addUndoState=False)
s.CoincidentConstraint(entity1=v[11], entity2=g[5], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.3136, 
    farPlane=4.20569, width=5.94658, height=2.96586, cameraPosition=(-3.952, 
    -4.17972, 3.75964), cameraTarget=(-3.952, -4.17972, 0))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-1']
s1 = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s1)
s2 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s2.geometry, s2.vertices, s2.dimensions, s2.constraints
s2.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s2, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.3579, 
    farPlane=4.16139, width=5.35602, height=2.67132, cameraPosition=(-3.35797, 
    -4.09395, 3.75964), cameraTarget=(-3.35797, -4.09395, 0))
s2.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__edit__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9232, 
    farPlane=28.6454, width=4.79256, height=2.39635, viewOffsetX=-4.09004, 
    viewOffsetY=-4.06451)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshEdgesInShaded=OFF, substructureTranslucency=OFF)
p = mdb.models['Model-1'].parts['Part-1']
s = p.features['Partition face-3'].sketch
mdb.models['Model-1'].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models['Model-1'].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Partition face-3'], filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.60799, 
    farPlane=3.9113, width=2.02184, height=1.00839, cameraPosition=(-4.81214, 
    -3.47053, 3.75964), cameraTarget=(-4.81214, -3.47053, 0))
s1.Line(point1=(-0.800845, 0.383046534805298), point2=(-0.800845, 0.52839))
s1.VerticalConstraint(entity=g[27], addUndoState=False)
s1.PerpendicularConstraint(entity1=g[2], entity2=g[27], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.35754, 
    farPlane=4.16175, width=3.8732, height=1.93176, cameraPosition=(-3.44201, 
    -4.04876, 3.75964), cameraTarget=(-3.44201, -4.04876, 0))
s1.radialPattern(geomList=(g[5], ), vertexList=(), number=33, totalAngle=-90.0, 
    centerPoint=(-0.800845, -0.800845))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.71103, 
    farPlane=3.80826, width=0.468292, height=0.233561, cameraPosition=(
    -4.90185, -3.70613, 3.75964), cameraTarget=(-4.90185, -3.70613, 0))
s1.Line(point1=(-0.800845, 0.52839), point2=(-0.73918836282329, 
    0.454204681916091))
s1.CoincidentConstraint(entity1=v[79], entity2=g[28], addUndoState=False)
s1.EqualDistanceConstraint(entity1=v[16], entity2=v[15], midpoint=v[79], 
    addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.72452, 
    farPlane=3.79477, width=0.468292, height=0.233561, cameraPosition=(
    -4.87996, -3.70316, 3.75964), cameraTarget=(-4.87996, -3.70316, 0))
s1.autoTrimCurve(curve1=g[60], point1=(-0.782930704803467, 0.506589320449829))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75151, 
    farPlane=3.76778, width=0.0783655, height=0.0390848, cameraPosition=(
    -4.9291, -3.73968, 3.75964), cameraTarget=(-4.9291, -3.73968, 0))
s1.autoTrimCurve(curve1=g[15], point1=(-0.741982313842773, 0.453866389541626))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.72979, 
    farPlane=3.7895, width=0.398048, height=0.198526, cameraPosition=(-4.8808, 
    -3.71759, 3.75964), cameraTarget=(-4.8808, -3.71759, 0))
s1.autoTrimCurve(curve1=g[61], point1=(-0.734487387390137, 0.454069283752442))
s1.autoTrimCurve(curve1=g[62], point1=(-0.664790960998535, 0.447163966445923))
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
p.features['Partition face-3'].setValues(sketch=s1)
del mdb.models['Model-1'].sketches['__edit__']
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
p = mdb.models['Model-1'].parts['Part-1']
p.regenerate()
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2126, 
    farPlane=28.3559, width=0.688768, height=0.344394, viewOffsetX=-4.74256, 
    viewOffsetY=-3.77249)
p = mdb.models['Model-1'].parts['Part-1']
del p.features['Partition face-3']
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.7894, 
    farPlane=28.7792, width=4.76959, height=2.38486, viewOffsetX=-3.81238, 
    viewOffsetY=-4.04963)
p = mdb.models['Model-1'].parts['Part-1']
f1, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f1[0], sketchPlaneSide=SIDE1, origin=(
    -4.199155, -4.199155, 0.0))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=3.75, 
    gridSpacing=0.09, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.74838, 
    farPlane=3.77091, width=0.108464, height=0.0540966, cameraPosition=(
    -4.99355, -3.75265, 3.75964), cameraTarget=(-4.99355, -3.75265, 0))
s.ArcByCenterEnds(center=(-0.800845, -0.800845), point1=(0.455718267402649, 
    -0.800845), point2=(-0.800845, 0.455718267402649), 
    direction=COUNTERCLOCKWISE)
s.CoincidentConstraint(entity1=v[10], entity2=g[3], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[1], entity2=v[3], midpoint=v[10], 
    addUndoState=False)
s.CoincidentConstraint(entity1=v[11], entity2=g[5], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[4], entity2=v[0], midpoint=v[11], 
    addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.68796, 
    farPlane=3.83133, width=0.690491, height=0.344383, cameraPosition=(
    -4.92946, -3.77863, 3.75964), cameraTarget=(-4.92946, -3.77863, 0))
s.Line(point1=(-0.800845, 0.52839), point2=(-0.800845, 0.383046534805298))
s.VerticalConstraint(entity=g[16], addUndoState=False)
s.PerpendicularConstraint(entity1=g[4], entity2=g[16], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.56957, 
    farPlane=3.94971, width=2.534, height=1.26383, cameraPosition=(-4.25198, 
    -3.97245, 3.75964), cameraTarget=(-4.25198, -3.97245, 0))
s.radialPattern(geomList=(g[5], ), vertexList=(), number=33, totalAngle=-90.0, 
    centerPoint=(-0.800845, -0.800845))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.6179, 
    farPlane=3.90138, width=1.36528, height=0.680935, cameraPosition=(-5.13069, 
    -3.61489, 3.75964), cameraTarget=(-5.13069, -3.61489, 0))
s.Line(point1=(-0.73918836282329, 0.454204681916091), point2=(-0.800845, 
    0.52839))
s.CoincidentConstraint(entity1=v[76], entity2=g[15], addUndoState=False)
#: Warning: Cannot continue yet--complete the step or cancel the procedure.
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.75318, 
    farPlane=3.76611, width=0.0622562, height=0.0310503, cameraPosition=(
    -4.82816, -3.75601, 3.75964), cameraTarget=(-4.82816, -3.75601, 0))
s.Line(point1=(-0.61646887559083, 0.442117867399294), point2=(
    -0.620027395935058, 0.438186791687012))
s.CoincidentConstraint(entity1=v[77], entity2=g[15], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.73948, 
    farPlane=3.77981, width=0.1942, height=0.0968574, cameraPosition=(-4.77388, 
    -3.75694, 3.75964), cameraTarget=(-4.77388, -3.75694, 0))
s.autoTrimCurve(curve1=g[50], point1=(-0.618300291748047, 0.43998923614502))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.64802, 
    farPlane=3.87126, width=1.48811, height=0.742193, cameraPosition=(-4.5874, 
    -3.74532, 3.75964), cameraTarget=(-4.5874, -3.74532, 0))
s.Line(point1=(-0.61646887559083, 0.442117867399294), point2=(
    -0.541524115814892, 0.502849122196788))
s.CoincidentConstraint(entity1=v[79], entity2=g[15], addUndoState=False)
s.Line(point1=(-0.541524115814892, 0.502849122196788), point2=(
    -0.495525031226671, 0.418060640996823))
s.CoincidentConstraint(entity1=v[80], entity2=g[15], addUndoState=False)
s.Line(point1=(-0.377521585066672, 0.382264686949949), point2=(
    -0.292168787780189, 0.42720801039764))
s.CoincidentConstraint(entity1=v[81], entity2=g[23], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[25], entity2=v[24], midpoint=v[81], 
    addUndoState=False)
s.Line(point1=(-0.292168787780189, 0.42720801039764), point2=(
    -0.263594974804599, 0.335074739864205))
s.CoincidentConstraint(entity1=v[82], entity2=g[15], addUndoState=False)
s.Line(point1=(-0.154842375975897, 0.276945264726674), point2=(
    -0.0623616013121887, 0.304373510108974))
s.CoincidentConstraint(entity1=v[83], entity2=g[27], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[33], entity2=v[32], midpoint=v[83], 
    addUndoState=False)
s.Line(point1=(-0.0623616013121887, 0.304373510108974), point2=(
    -0.0523111355575023, 0.208436080159735))
s.CoincidentConstraint(entity1=v[84], entity2=g[15], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.60515, 
    farPlane=3.91414, width=1.48811, height=0.742194, cameraPosition=(-4.15064, 
    -4.05667, 3.75964), cameraTarget=(-4.15064, -4.05667, 0))
s.Line(point1=(0.0430113145560775, 0.130206967061697), point2=(
    0.139066082290501, 0.139066082290501))
s.CoincidentConstraint(entity1=v[85], entity2=g[15], addUndoState=False)
s.Line(point1=(0.139066082290501, 0.139066082290501), point2=(
    0.130206967061697, 0.0430113145560772))
s.CoincidentConstraint(entity1=v[86], entity2=g[15], addUndoState=False)
s.Line(point1=(0.208436080159735, -0.0523111355575026), point2=(
    0.304373510108974, -0.062361601312189))
s.CoincidentConstraint(entity1=v[87], entity2=g[15], addUndoState=False)
s.Line(point1=(0.304373510108974, -0.062361601312189), point2=(
    0.276945264726675, -0.154842375975897))
s.CoincidentConstraint(entity1=v[88], entity2=g[37], addUndoState=False)
s.EqualDistanceConstraint(entity1=v[53], entity2=v[52], midpoint=v[88], 
    addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58714, 
    farPlane=3.93215, width=1.66161, height=0.828727, cameraPosition=(-3.96265, 
    -4.34044, 3.75964), cameraTarget=(-3.96265, -4.34044, 0))
s.Line(point1=(0.335074739864205, -0.2635949748046), point2=(0.427208010397641, 
    -0.29216878778019))
s.CoincidentConstraint(entity1=v[89], entity2=g[15], addUndoState=False)
s.Line(point1=(0.427208010397641, -0.29216878778019), point2=(0.38226468694995, 
    -0.377521585066673))
s.CoincidentConstraint(entity1=v[90], entity2=g[15], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.59592, 
    farPlane=3.92337, width=1.57704, height=0.786547, cameraPosition=(-3.8958, 
    -4.49838, 3.75964), cameraTarget=(-3.8958, -4.49838, 0))
s.Line(point1=(0.418060640996822, -0.495525031226671), point2=(
    0.502849122196788, -0.541524115814892))
s.CoincidentConstraint(entity1=v[91], entity2=g[15], addUndoState=False)
s.Line(point1=(0.502849122196788, -0.541524115814892), point2=(
    0.442117867399294, -0.61646887559083))
s.CoincidentConstraint(entity1=v[92], entity2=g[15], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(width=1.34048, height=0.668565, 
    cameraPosition=(-3.8188, -4.71228, 3.75964), cameraTarget=(-3.8188, 
    -4.71228, 0))
s.Line(point1=(0.454204681916091, -0.739188362823291), point2=(0.52839, 
    -0.800845))
s.CoincidentConstraint(entity1=v[93], entity2=g[15], addUndoState=False)
s.CoincidentConstraint(entity1=v[75], entity2=g[3], addUndoState=False)
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.53304, 
    farPlane=3.98625, width=2.18275, height=1.08865, cameraPosition=(-4.88432, 
    -3.96732, 3.75964), cameraTarget=(-4.88432, -3.96732, 0))
s.autoTrimCurve(curve1=g[17], point1=(-0.743116709442138, 0.492163327484131))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.6591, 
    farPlane=3.86019, width=0.968497, height=0.483038, cameraPosition=(
    -4.92633, -3.84391, 3.75964), cameraTarget=(-4.92633, -3.84391, 0))
s.autoTrimCurve(curve1=g[15], point1=(-0.762070509643554, 0.455150035171509))
s.autoTrimCurve(curve1=g[19], point1=(-0.609633299560547, 0.468351033477783))
s.autoTrimCurve(curve1=g[21], point1=(-0.488417956085205, 0.452149775772095))
s.autoTrimCurve(curve1=g[67], point1=(-0.520252081604004, 0.426947739868164))
s.autoTrimCurve(curve1=g[71], point1=(-0.567391726226806, 0.43414845779419))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.64135, 
    farPlane=3.87793, width=1.57704, height=0.786547, cameraPosition=(-4.98231, 
    -3.90857, 3.75964), cameraTarget=(-4.98231, -3.90857, 0))
s.autoTrimCurve(curve1=g[23], point1=(-0.362474772186279, 0.41906991317749))
s.autoTrimCurve(curve1=g[25], point1=(-0.245841833801269, 0.37412419631958))
s.autoTrimCurve(curve1=g[70], point1=(-0.289703699798584, 0.347743180541992))
s.autoTrimCurve(curve1=g[76], point1=(-0.342537733764648, 0.368261960296631))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.65587, 
    farPlane=3.86342, width=0.999627, height=0.498564, cameraPosition=(
    -4.46222, -4.00396, 3.75964), cameraTarget=(-4.46222, -4.00396, 0))
s.autoTrimCurve(curve1=g[29], point1=(-0.0305322138977049, 0.233280089645386))
s.autoTrimCurve(curve1=g[27], point1=(-0.135423514099121, 0.305122998504639))
s.autoTrimCurve(curve1=g[75], point1=(-0.120890471191406, 0.256195452957154))
s.autoTrimCurve(curve1=g[80], point1=(-0.0855053393554686, 0.235138085632324))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.63755, 
    farPlane=3.88174, width=1.62773, height=0.811829, cameraPosition=(-4.55486, 
    -4.02895, 3.75964), cameraTarget=(-4.55486, -4.02895, 0))
s.autoTrimCurve(curve1=g[31], point1=(0.068201211242676, 0.153560307769776))
s.autoTrimCurve(curve1=g[33], point1=(0.15154233291626, 0.06279721572876))
s.autoTrimCurve(curve1=g[82], point1=(0.105241921691895, 0.0638057263183596))
s.autoTrimCurve(curve1=g[86], point1=(0.071288255004883, 0.0970856221008303))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.63755, 
    farPlane=3.88174, width=1.17603, height=0.586546, cameraPosition=(-4.08982, 
    -4.30992, 3.75964), cameraTarget=(-4.08982, -4.30992, 0))
s.autoTrimCurve(curve1=g[37], point1=(0.301127818374634, -0.142464968414306))
s.autoTrimCurve(curve1=g[35], point1=(0.226046231536865, -0.038999888153076))
s.autoTrimCurve(curve1=g[85], point1=(0.230506566314697, -0.084903094024658))
s.autoTrimCurve(curve1=g[90], point1=(0.259498503952027, -0.120606276245117))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58523, 
    farPlane=3.93406, width=1.68005, height=0.837923, cameraPosition=(-4.04722, 
    -4.39771, 3.75964), cameraTarget=(-4.04722, -4.39771, 0))
s.autoTrimCurve(curve1=g[39], point1=(0.370699075012207, -0.249042841644287))
s.autoTrimCurve(curve1=g[41], point1=(0.43760338142395, -0.362500521392822))
s.autoTrimCurve(curve1=g[92], point1=(0.353707459716797, -0.30316958114624))
s.autoTrimCurve(curve1=g[95], point1=(0.374946740417481, -0.351050707550049))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.58523, 
    farPlane=3.93406, height=0.837923, cameraPosition=(-3.93205, -4.63332, 
    3.75964), cameraTarget=(-3.93205, -4.63332, 0))
s.autoTrimCurve(curve1=g[43], point1=(0.457194951324463, -0.490895601959228))
s.autoTrimCurve(curve1=g[45], point1=(0.494363930969238, -0.607476088256836))
s.autoTrimCurve(curve1=g[97], point1=(0.435955193786621, -0.590822073669433))
s.autoTrimCurve(curve1=g[101], point1=(0.426397469787598, -0.532531592102051))
s.autoTrimCurve(curve1=g[47], point1=(0.489054110794068, -0.73758873626709))
s.autoTrimCurve(curve1=g[100], point1=(0.455070880203247, -0.772979113311767))
session.viewports['Viewport: 1'].view.setValues(nearPlane=3.47563, 
    farPlane=4.04365, width=3.7864, height=1.88847, cameraPosition=(-4.15041, 
    -4.37427, 3.75964), cameraTarget=(-4.15041, -4.37427, 0))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#1 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#f8000000 #ffff ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #e000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0986, 
    farPlane=28.47, width=2.46711, height=1.23359, viewOffsetX=-4.33067, 
    viewOffsetY=-3.89254)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1c00 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.8603, 
    farPlane=28.7083, width=5.6256, height=2.81288, viewOffsetX=-3.70324, 
    viewOffsetY=-3.72053)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #20000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshEdgesInShaded=ON, substructureTranslucency=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1231, 
    farPlane=28.4455, width=2.14171, height=1.07088, viewOffsetX=-4.40904, 
    viewOffsetY=-3.78229)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #200 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #80 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #100 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #40 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #20 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #2 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #4 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#80000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#40000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#20000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#8000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#10000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.9779, 
    farPlane=28.5906, width=2.949, height=1.47454, viewOffsetX=-3.84567, 
    viewOffsetY=-3.9179)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#7f80000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#7f800 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0068, 
    farPlane=28.5617, width=3.07377, height=1.53693, viewOffsetX=-3.70286, 
    viewOffsetY=-4.01805)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#7e0 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0009, 
    farPlane=28.5677, width=3.13981, height=1.56995, viewOffsetX=-3.48917, 
    viewOffsetY=-4.21957)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1f ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2195, 
    farPlane=28.3491, width=0.622974, height=0.311496, viewOffsetX=-4.61998, 
    viewOffsetY=-3.77108)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1400 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:3 #1000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #1000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #400 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #600 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:3 #100 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #400 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #200 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.2274, 
    farPlane=28.3412, width=0.546702, height=0.273359, viewOffsetX=-4.41002, 
    viewOffsetY=-3.85241)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #18 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:2 #2000000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1426, 
    farPlane=28.426, width=1.36265, height=0.681344, viewOffsetX=-3.94039, 
    viewOffsetY=-3.936)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#60000000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:2 #400 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#40000000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#20000000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.0547, 
    farPlane=28.5138, width=2.64695, height=1.32352, viewOffsetX=-3.30245, 
    viewOffsetY=-4.48173)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#60 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#4000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#40 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#20 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#14 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#80 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#10 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1402, 
    farPlane=28.4284, width=1.38593, height=0.692987, viewOffsetX=-3.5291, 
    viewOffsetY=-4.64878)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#30 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#1000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#1 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#35 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=28.1149, 
    farPlane=28.4536, width=1.62905, height=0.814547, viewOffsetX=-3.60194, 
    viewOffsetY=-4.37596)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#60000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0 #1000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=1, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#60000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=24.4777, 
    farPlane=32.0908, width=37.3391, height=18.6701, viewOffsetX=2.73162, 
    viewOffsetY=-5.25077)
elemType1 = mesh.ElemType(elemCode=CPE8, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=CPE6M, elemLibrary=STANDARD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.0959, 
    farPlane=30.4727, width=28.7609, height=14.3808, viewOffsetX=-1.90838, 
    viewOffsetY=-4.8782)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
session.viewports['Viewport: 1'].view.setValues(nearPlane=20.0407, 
    farPlane=36.5278, width=98.4268, height=49.2148, viewOffsetX=0.123212, 
    viewOffsetY=-11.3284)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.3866, 
    farPlane=30.1819, width=24.9793, height=12.49, viewOffsetX=2.91817, 
    viewOffsetY=-0.906829)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:3 #12000000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=10, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#0:3 #8000000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=8, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges1 = e.getSequenceFromMask(mask=('[#0:3 #10000000 ]', ), )
pickedEdges2 = e.getSequenceFromMask(mask=('[#0:3 #2000000 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, 
    end2Edges=pickedEdges2, ratio=10.0, number=10, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.4593, 
    farPlane=30.1093, width=24.0322, height=12.0165, viewOffsetX=2.22409, 
    viewOffsetY=-1.60358)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.setMeshControls(regions=pickedRegions, elemShape=QUAD)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.generateMesh(regions=pickedRegions)
#: 
#: Part: Part-1
#:   Number of elements :  721,   Analysis errors:  0 (0%),  Analysis warnings:  43 (5.96394%)
#: 
#: Part: Part-1
#:   Number of elements :  721,   Analysis errors:  0 (0%),  Analysis warnings:  43 (5.96394%)
#: 
#: Part: Part-1
#:   Number of elements :  721,   Analysis errors:  0 (0%),  Analysis warnings:  43 (5.96394%)
session.viewports['Viewport: 1'].view.setValues(nearPlane=25.5683, 
    farPlane=31.0002, width=35.5909, height=17.796, viewOffsetX=8.21524, 
    viewOffsetY=-3.46844)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.deleteMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=25.0777, 
    farPlane=31.4908, width=41.9063, height=20.9538, viewOffsetX=7.87557, 
    viewOffsetY=-1.8028)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=FREE)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=25.9022, 
    farPlane=30.6663, width=31.2728, height=15.6368, viewOffsetX=4.63327, 
    viewOffsetY=0.966774)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.setMeshControls(regions=pickedRegions, elemShape=QUAD_DOMINATED, 
    technique=STRUCTURED)
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.3433, 
    farPlane=29.2253, width=12.4507, height=6.22553, viewOffsetX=-1.0225, 
    viewOffsetY=-2.48169)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.3868, 
    farPlane=30.1817, width=24.977, height=12.4888, viewOffsetX=1.57561, 
    viewOffsetY=-1.3661)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#8a28a28a #150a2 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f, e, d1 = p.faces, p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=f[48], sketchPlaneSide=SIDE1, origin=(
    0.062422, 0.062422, 0.0))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=28.28, gridSpacing=0.7, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts['Part-1']
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.Line(point1=(-5.062422, -5.062422), point2=(4.937578, 4.937578))
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.6696, 
    farPlane=28.899, width=8.1952, height=4.08736, cameraPosition=(-2.0435, 
    -3.6463, 28.2843), cameraTarget=(-2.0435, -3.6463, 0))
s1.autoTrimCurve(curve1=g[249], point1=(-4.35416575534058, -4.38148900134277))
s1.autoTrimCurve(curve1=g[250], point1=(-4.97579928500366, -4.98570748431396))
s1.autoTrimCurve(curve1=g[251], point1=(-4.15213509661865, -4.15300341708374))
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedFaces = f.getSequenceFromMask(mask=('[#0 #10000 ]', ), )
e1, d2 = p.edges, p.datums
p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
session.viewports['Viewport: 1'].view.setValues(nearPlane=27.5481, 
    farPlane=29.0205, width=7.10552, height=3.55286, viewOffsetX=-2.3696, 
    viewOffsetY=-3.32432)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #8000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#0 #2540 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#10000000 #5 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#4514514 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.0746, 
    farPlane=30.4939, width=29.037, height=14.5189, viewOffsetX=2.31999, 
    viewOffsetY=-1.33872)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges1 = e.getSequenceFromMask(mask=('[#1 ]', ), )
p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges1, ratio=10.0, 
    number=10, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 #20000 ]', ), )
p.generateMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 #20000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 #20000 ]', ), )
p.setMeshControls(regions=pickedRegions, technique=STRUCTURED)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 #20000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=25.5199, 
    farPlane=31.0487, width=36.2159, height=18.1085, viewOffsetX=-1.17728, 
    viewOffsetY=-1.50013)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 #20000 ]', ), )
p.deleteMesh(regions=pickedRegions)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
pickedEdges = e.getSequenceFromMask(mask=('[#400 #0:2 #10000000 ]', ), )
p.seedEdgeByNumber(edges=pickedEdges, number=8, constraint=FINER)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
pickedRegions = f.getSequenceFromMask(mask=('[#1 #20000 ]', ), )
p.generateMesh(regions=pickedRegions)
session.viewports['Viewport: 1'].view.setValues(nearPlane=25.2677, 
    farPlane=31.3008, width=39.4646, height=19.7329, viewOffsetX=1.44659, 
    viewOffsetY=-0.591037)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON, mesh=OFF)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((210000.0, 0.3), 
    ))
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
    material='Material-1', thickness=None)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, 
    adaptiveMeshConstraints=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.Job(name='Mesh_720', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
mdb.jobs['Mesh_720'].writeInput(consistencyChecking=OFF)
#: The job input file has been written to "Mesh_720.inp".
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=26.7772, 
    farPlane=29.7914, width=19.8805, height=9.94052, viewOffsetX=0.276552, 
    viewOffsetY=-1.49838)
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
mdb.save()
#: The model database has been saved to "/disk501/home2/dp44dafy/PPP/verification_model/cae_model/mesh_gen.cae".
