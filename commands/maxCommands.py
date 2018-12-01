from Py3dsMax import mxs
import numpy as np
from contextlib import contextmanager

def thingIsMesh(thing):
	isMesh = False
	if mxs.classOf(thing) == mxs.XRefObject:
		thing = thing.actualBaseObject
	if mxs.classOf(thing) == mxs.Editable_Mesh:
		isMesh = True
	return isMesh

@contextmanager
def polyManager(thing):
	"""
		So max "knows" what a polygon is, but for meshes
		everithing is a triangle. We need to work in polygons for
		cross ddc compatibility.
		Because mxs.getPolysUsingFace() is *SLOW*, we just add
		and remove a poly select using this context manager
	"""
	isMesh = thingIsMesh(thing)
	selectModifier = None
	mxs.disableSceneRedraw()
	if isMesh:
		selectModifier = mxs.poly_select()
		mxs.addmodifier(thing, selectModifier)
	try:
		yield thing
	finally:
		if isMesh:
			mxs.deleteModifier(thing, selectModifier)
	mxs.enableSceneRedraw()
	mxs.redrawViews()



def getSingleSelection():
	return mxs.selection[0]

def getObjectByName(name):
	return mxs.getNodeByName(name)

def getObjectName(thing):
	return thing.name

def getFaces(thing):
	topo = []
	with polyManager(thing) as obj:
		for i in xrange(mxs.getNumFaces(obj)):
			fv = mxs.polyop.getFaceVerts(obj, i+1)
			topo.append([x-1 for x in fv])
	return topo

def getVerts(thing):
	isMesh = thingIsMesh(thing)
	op = mxs if isMesh else mxs.polyop
	maxVerts = [op.getVert(thing, i + 1) for i in xrange(op.getNumVerts(thing))]
	return np.array([(p.x, p.y, p.z) for p in maxVerts])

def createRawObject(name, faces, verts, uvFaces, uvs):
	mm = mxs.mesh(numverts=len(verts))
	oldUndo = mxs.execute("set undo off")
	try:
		mm.name = name
		for i in xrange(len(verts)):
			mxs.meshop.setvert(mm, i+1, mxs.Point3(*verts[i]))

		mxs.convertTo(mm, mxs.PolyMeshObject)
		for face in faces:
			fPlus = [i+1 for i in face]
			mxs.polyop.createPolygon(mm, fPlus)

		if uvs is not None:
			mxs.polyop.defaultMapFaces(mm, 1)
			mxs.polyop.setNumMapVerts(mm, 1, len(uvs), keep=False)
			mxs.polyop.setNumMapFaces(mm, 1, len(uvFaces), keep=False)

			for i in xrange(len(uvs)):
				mxs.polyop.setMapVert(mm, 1, i+1, mxs.Point3(*uvs[i]))

			mxs.convertTo(mm, mxs.PolyMeshObject)
			for f, face in enumerate(uvFaces):
				fPlus = [i+1 for i in face]
				mxs.polyop.setMapFace(mm, 1, f+1, fPlus)
	finally:
		mxs.execute("set undo {}".format('on' if oldUndo else 'off'))
	return mm

def selectVerts(obj, idx):
	selMod = None
	if list(obj.modifiers):
		if obj.modifiers[0] == mxs.Mesh_Select:
			selMod = obj.modifiers[0]
		elif obj.modifiers[0] == mxs.Poly_Select:
			selMod = obj.modifiers[0]

	if selMod is None:
		if mxs.classOf(obj) == mxs.Editable_Mesh:
			selMod = mxs.addModifier(obj, mxs.Mesh_Select)
		elif mxs.classOf(obj) == mxs.Editable_Poly:
			selMod = mxs.addModifier(obj, mxs.Poly_Select)

	if selMod is None:
		raise RuntimeError("Could not build selection")

	mxs.select(obj)
	mxs.subObjectLevel = 1
	sel = [] if idx is None else [idx + 1]
	obj.selectedVerts = sel
	mxs.completeRedraw()

def getVertSelection(obj):
	return obj.selectedVerts[0]

