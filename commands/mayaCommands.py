import re
import maya.OpenMaya as om
import maya.cmds as cmds
import numpy as np
from ctypes import c_float
from itertools import groupby, count
from operator import itemgetter

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

def getUVs(thing):
	# Get the MDagPath from the name of the mesh
	sl = om.MSelectionList()
	sl.add(thing)
	thing = om.MDagPath()
	sl.getDagPath(0, thing)
	meshFn = om.MFnMesh(thing)

	vIdx = om.MIntArray()
	util = om.MScriptUtil()
	util.createFromInt(0)
	uvIdxPtr = util.asIntPtr()
	uArray = om.MFloatArray()
	vArray = om.MFloatArray()
	meshFn.getUVs(uArray, vArray)
	hasUVs = uArray.length() > 0
	if not hasUVs:
		return None, None

	uvFaces = []
	for i in range(meshFn.numPolygons()):
		meshFn.getPolygonVertices(i, vIdx)
		uvFace = []
		for j in reversed(xrange(vIdx.length())):
			meshFn.getPolygonUVid(i, j, uvIdxPtr)
			uvIdx = util.getInt(uvIdxPtr)
			if uvIdx >= uArray.length() or uvIdx < 0:
				uvIdx = 0
			uvFace.append(uvIdx)
		uvFaces.append(uvFace)

	uvs = [(uArray[i], vArray[i]) for i in xrange(uArray.length())]
	return uvs, uvFaces

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
		vertArray.set(i, float(v[0]), float(v[1]), float(v[2]))

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
	sel = cmds.ls(selection=True, long=True)
	sel = [i for i in sel if '.vtx' in i]
	longName = cmds.ls(obj, long=True)[0]
	nums = []
	for s in sel:
		if s.startswith(longName):
			nums.append(int(s.split('[')[-1][:-1]))
	return nums

def cloneObject(obj, name):
	return cmds.duplicate(obj, name=name)[0]

def freezeObject(obj):
	cmds.delete(obj, constructionHistory=True)

def setObjectName(obj, newName):
	cmds.rename(obj, newName)

def setAllVerts(obj, newVerts):
	# Get the api objects
	sel = om.MSelectionList()
	sel.add(obj)
	dagPath = om.MDagPath()
	sel.getDagPath(0, dagPath)
	fnMesh = om.MFnMesh(dagPath)

	# Build the scriptUtil object
	ptCount = fnMesh.numVertices()
	util = om.MScriptUtil()
	util.createFromList([0.0] * (4 * ptCount), 4 * ptCount)
	ptr = util.asFloat4Ptr()

	# Get the scriptUtil object as a pointer in numpy
	cta = (c_float * 4 * ptCount).from_address(int(ptr))
	out = np.ctypeslib.as_array(cta)

	# Copy the input values into the pointer
	out[:, :3] = newVerts

	# Set the vert positions
	fnMesh.setPoints(om.MFloatPointArray(ptr, ptCount))

def selectAdjacentEdges(obj, centers):
	sel = []
	for k, g in groupby(enumerate(centers), lambda (i,x):i-x):
		adj = list(map(itemgetter(1), g))
		first, last = adj[0], adj[-1]
		if first == last:
			sel.append('{0}.vtx[{1}]'.format(obj, first))
		else:
			sel.append('{0}.vtx[{1}:{2}]'.format(obj, first, last))

	edges = cmds.polyListComponentConversion(*sel, fromVertex=True, toEdge=True)
	cmds.select(*edges, r=True)

def rootWindow():
	"""
	Returns the currently active QT main window
	Only works for QT UI's like Maya
	"""
	from MeshCrawler.Qt.QtWidgets import QApplication, QSplashScreen, QDialog, QMainWindow
	# for MFC apps there should be no root window
	window = None
	if QApplication.instance():
		inst = QApplication.instance()
		window = inst.activeWindow()
		# Ignore QSplashScreen's, they should never be considered the root window.
		if isinstance(window, QSplashScreen):
			return None
		# If the application does not have focus try to find A top level widget
		# that doesn't have a parent and is a QMainWindow or QDialog
		if window is None:
			windows = []
			dialogs = []
			for w in QApplication.instance().topLevelWidgets():
				if w.parent() is None:
					if isinstance(w, QMainWindow):
						windows.append(w)
					elif isinstance(w, QDialog):
						dialogs.append(w)
			if windows:
				window = windows[0]
			elif dialogs:
				window = dialogs[0]

		# grab the root window
		if window:
			while True:
				parent = window.parent()
				if not parent:
					break
				if isinstance(parent, QSplashScreen):
					break
				window = parent

	return window

