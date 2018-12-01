import re
import maya.OpenMaya as om
import maya.cmds as cmds
import numpy as np
from ctypes import c_float

def getSingleSelection():
	thing = cmds.ls(sl=True, objectsOnly=True)
	if thing:
		return thing[0]
	return None

def getObjectByName(name):
	thing = cmds.ls(name, objectsOnly=True)
	if thing:
		return thing[0]
	return None

def getObjectName(thing):
	return thing

def getFaces(thing):
	sel = om.MSelectionList()
	sel.add(thing)
	dagPath = om.MDagPath()
	sel.getDagPath(0, dagPath)
	fnMesh = om.MFnMesh(dagPath)
	faces = []
	vIdx = om.MIntArray()
	for i in range(fnMesh.numPolygons()):
		fnMesh.getPolygonVertices(i, vIdx)
		face = [vIdx[j] for j in xrange(vIdx.length())]
		faces.append(face)
	return faces

def getVerts(thing):
	sel = om.MSelectionList()
	sel.add(thing)
	dagPath = om.MDagPath()
	sel.getDagPath(0, dagPath)

	fnMesh = om.MFnMesh(dagPath)
	# rawPts is a SWIG float pointer
	rawPts = fnMesh.getRawPoints()
	ptCount = fnMesh.numVertices()

	cta = (c_float * 3 * ptCount).from_address(int(rawPts))
	out = np.ctypeslib.as_array(cta)

	# for safety, make a copy of out so I don't corrupt memory
	return np.copy(out)

def createRawObject(name, faces, verts, uvFaces, uvs):
	dup = cmds.polyPlane(name=name, constructionHistory=False)[0]
	sel = om.MSelectionList()
	sel.add(dup)
	dagPath = om.MDagPath()
	sel.getDagPath(0, dagPath)
	fnMesh = om.MFnMesh(dagPath)

	counts = om.MIntArray()
	connects = om.MIntArray()
	for face in faces:
		counts.append(len(face))
		for f in face:
			connects.append(f)

	vertArray = om.MFloatPointArray()
	vertArray.setLength(len(verts))
	for i, v in enumerate(verts):
		vertArray.set(i, v[0], v[1], v[2])

	fnMesh.createInPlace(len(verts), len(faces), vertArray, counts, connects)
	fnMesh.updateSurface()

	if uvs is not None:
		us = om.MFloatArray()
		vs = om.MFloatArray()
		for uv in uvs:
			us.append(uv[0])
			vs.append(uv[1])
		fnMesh.setUvs(us, vs)

		uvConnects = om.MIntArray()
		for uvFace in uvFaces:
			for uvf in uvFace:
				uvConnects.append(uvf)

		fnMesh.assignUVs(counts, uvConnects)
	cmds.sets(dup, e=True, fe='initialShadingGroup')
	return dup

def selectVerts(obj, idx):
	cmds.selectType(objectComponent=True, allComponents=False)
	cmds.selectType(objectComponent=True, vertex=True)
	cmds.selectType(vertex=True)
	cmds.hilite(obj)
	if idx is not None:
		cmds.select('pSphere1.vtx[{0}]'.format(idx), replace=True)

def getVertSelection(obj):
	sel = cmds.ls(selection=True)
	sel = [i for i in sel if '.vtx' in i]
	longs = cmds.ls(sel, long=True)
	for s in longs:
		if s.startswith(obj):
			return int(s.split('[')[:-1])
	return None
