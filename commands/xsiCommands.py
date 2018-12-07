import numpy as np
from dcc.xsi import xsi
from itertools import chain

def getSingleSelection():
	return xsi.Selection[0]

def getObjectByName(name):
	return xsi.GetValue(name)

def getObjectName(thing):
	return thing.Name

def getFaces(thing):
	vertArray, faceArray = thing.ActivePrimitive.Geometry.Get2()
	ptr = 0
	faces = []
	while ptr < len(faceArray):
		count = faceArray[ptr]
		ptr += 1
		indices = list(reversed(faceArray[ptr:ptr+count]))
		ptr += count
		faces.append(indices)
	return faces

def getVerts(thing):
	vertArray, faceArray = thing.ActivePrimitive.Geometry.Get2()
	verts = np.array(vertArray)
	return verts.T

def _getTexProps(thing, name=None):
	textureCls = [cluster for cluster in thing.ActivePrimitive.Geometry.Clusters if cluster.Type == "sample"]
	out = []
	for cluster in textureCls:
		spaces = [i for i in cluster.Properties if i.Type == 'uvspace']
		if spaces:
			out.extend(spaces)
	if name is not None:
		out = [i for i in out if i.name == name]
	return out

def _getUVs(thing):
	""" Get the direct xsi uvws (without indexing) """
	texProps = _getTexProps(thing)
	if not texProps:
		return None
	texProp = texProps[0]
	return  zip(*texProp.Elements.Array)

def getUVs(thing):
	bigUvs = _getUVs(thing)
	if bigUvs is None:
		return None, None

	# Make a dict that keeps track of the first index
	# each uv value is found at
	uvd = {}
	count = 0
	for uv in bigUvs:
		if uv not in uvd:
			uvd[uv] = count
			count += 1

	# Build the collapsed uv list
	uvs = [None] * len(uvd)
	for uv, idx in uvd.iteritems():
		uvs[idx] = uv
	uvs = np.array(uvs)

	# loop over the faces
	faces = getFaces(thing)
	ptr = 0
	uvFaces = []
	for face in faces:
		uvFace = [uvd[fuv] for fuv in reversed(bigUvs[ptr: ptr+len(face)])]
		uvFaces.append(uvFace)
		ptr += len(face)

	return uvs, uvFaces

def _createRawObject(name, faces, verts, uvws=None):
	''' After converting the values to xsi style, this does the actual creation '''
	dup = xsi.ActiveSceneRoot.AddPolygonMesh(verts, faces, name)

	if uvws is not None:
		texName = "Texture_Projection"
		xsi.CreateProjection(dup, "", "", "", texName, True, "", "")

		texProps = _getTexProps(dup, name=texName)
		if texProps:
			texProp = texProps[0]
			xsi.FreezeObj(texProp)
			if len(uvws[0]) == texProp.Elements.Count:
				texProp.Elements.Array = uvws
			xsi.FreezeObj(texProp)

	return dup

def createRawObject(name, faces, verts, uvFaces, uvs):
	vertArray = verts.T.tolist()
	faceArray = []
	for face in faces:
		faceArray.append(len(face))
		faceArray.extend(reversed(face))

	uvws = None
	if uvs is not None:
		uvws = [(i[0], i[1], 0.0) for i in uvs]
		uvws = [uvws[i] for i in chain.from_iterable(map(reversed, uvFaces))]
		uvws = zip(*uvws)

	return _createRawObject(name, faceArray, vertArray, uvws)

def selectVerts(obj, idx):
	subComp = obj.ActivePrimitive.Geometry.Points.SubComponent
	xsi.SelectFilter("Object")
	if idx is None:
		xsi.DeselectAll()
		return

	subComp.ElementArray = [idx]
	cc = subComp.ComponentCollection
	xsi.SelectGeometryComponents(cc)

def getVertSelection(obj):
	sel = xsi.selection
	for i in sel:
		if i.type != 'pntSubComponent':
			continue
		sc = i.SubComponent
		thing = sc.Parent3DObject
		if thing.IsEqualTo(obj):
			return list(sc.ElementArray)
	return None

def cloneObject(obj, name):
	vertArray, faceArray = obj.ActivePrimitive.Geometry.Get2()
	uvws = _getUVs(obj)
	if uvws is not None:
		uvws = zip(*uvws)
	return _createRawObject(name, faceArray, vertArray, uvws=uvws)

def freezeObject(obj):
	xsi.FreezeObj(obj)

def setObjectName(obj, newName):
	obj.Name = newName

def setAllVerts(obj, newVerts):
	verts = newVerts.T.tolist()
	obj.ActivePrimitive.Geometry.Points.PositionArray = verts

def selectAdjacentEdges(obj, centers):
	subComp = obj.ActivePrimitive.Geometry.Points.SubComponent
	subComp.ElementArray = sorted(list(centers))
	cc = subComp.ComponentCollection
	edges = [edge.Index for edge in cc.NeighborEdges()]
	subComp = obj.ActivePrimitive.Geometry.Edges.SubComponent
	subComp.ElementArray = edges
	xsi.SelectGeometryComponents(subComp)

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

