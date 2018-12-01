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
		indices = reversed(faceArray[ptr:ptr+count])
		ptr += count
		faces.append(indices)
	return faces

def getVerts(thing):
	vertArray, faceArray = thing.ActivePrimitive.Geometry.Get2()
	verts = np.array(vertArray)
	return verts.T

def createRawObject(name, faces, verts, uvFaces, uvs):
	vertArray = verts.T.tolist()

	faceArray = []
	for face in faces:
		faceArray.append(len(face))
		faceArray.extend(reversed(face))

	dup = xsi.ActiveSceneRoot.AddPolygonMesh(vertArray, faceArray, name)

	if uvs is not None:
		uvws = [(i[0], i[1], 0.0) for i in uvs]
		uvws = [uvws[i] for i in chain.from_iterable(map(reversed, uvFaces))]

		texName = "Texture_Projection"
		xsi.CreateProjection(dup, "", "", "", texName, True, "", "")

		texProp = None
		textureCls = [cluster for cluster in dup.ActivePrimitive.Geometry.Clusters if cluster.Type == "sample"]
		for cluster in textureCls:
			texProp = cluster.Properties(texName)
			if texProp:
				break

		if texProp is not None:
			xsi.FreezeObj(texProp)
			if len(uvws) == texProp.Elements.Count:
				texProp.Elements.Array = zip(*uvws)
			xsi.FreezeObj(texProp)

	return dup

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
			ea = list(sc.ElementArray)
			if ea:
				return ea[0]
	return None

