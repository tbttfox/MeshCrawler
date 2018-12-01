#pylint:disable=missing-docstring
'''
Meshcrawler is an awesome library that allows you to build a 1:1 topology match
between two meshes that don't have the same vertex order.

This is done (in most cases) *without* users specifying any matching vertex points.
It uses some relatively simple heuristics to accomplish this, and they exist mostly
in meshcrawlergen.py

This is the actual matching algorithm. It takes meshes (that provide adjacency
info), and matched ordered pairs, and either returns a set of vertex matches,
or raises an error
'''

import blurdev
import sys, copy, time
from meshcrawlerErrors import TopologyMismatch

from Qt.QtWidgets import QApplication

###################################################
###				  Match Point Order				###
###################################################

def flipMultiDict(e, f):
	""" Make a dict keyed off the values of two dicts

	Args:
		e,f: Two dictionaries with hashable values

	Returns:
		A dictionary in this form: (eValueSet,fValueSet):key
	"""
	inv = {}
	allkeys = set(e.keys()) | set(f.keys())
	for k in allkeys:
		eVal = frozenset(e.get(k, ()))
		fVal = frozenset(f.get(k, ()))
		inv.setdefault((eVal, fVal), []).append(k)
	return inv

def growTracked(mesh, growSet, allSet):
	""" Grow a set of verts along edges and faces

	This function takes a set of vertices and returns two
	dicts. One that is vert:(edgeAdjVerts..), and another
	that is vert:(faceAdjVerts..)

	While I'm growing a vert, if all vertices adjacent to
	that vert are matched, then I no longer need to grow
	from there.

	Args:
		growSet: A set of Vertex objects to grow. This set
			will have all useless verts removed and will be
			changed in-place.
		allSet: A set of Vertex objects to exclude from
			the growth

	Returns:
		edgeDict: A dictionary with verts as keys and all
			edge adjacent vertices as the value
		faceDict: A dictionary with verts as keys and all
			face adjacent vertices as the value
		newGrowSet: The original growSet with all used-up
			vertices removed
	"""
	edgeDict = {}
	faceDict = {}
	# Verts that are grown, but have no adjacent verts
	# that aren't in the allSet
	rem = []
	for vert in growSet:
		edgeFound = False
		faceFound = False
		for eadj in mesh.adjacentVertsByEdge(vert):
			if eadj not in allSet:
				edgeFound = True
				edgeDict.setdefault(eadj, []).append(vert)

		for fadj in mesh.adjacentVertsByFace(vert):
			if fadj not in allSet:
				faceFound = True
				faceDict.setdefault(fadj, []).append(vert)

		if not edgeFound and not faceFound:
			rem.append(vert)

	newGrowSet = growSet - set(rem)
	return edgeDict, faceDict, newGrowSet

def matchByTopology(orderMesh, shapeMesh, vertexPairs, matchedNum=None, vertNum=None, symmetry=False, pBar=None):
	""" Match the topology of two meshes with different vert orders

	Provide a 1:1 vertex index match between two meshes that don't
	necessarily have the same vertex order.
	At minimum, 3 vertex pairs around a single poly are required.

	The algorithm simultaneously performs two different types of
	"grow verts" operations on each mesh and each vertex keeps track
	of where it was grown from, and how.
	New matching verts will be grown from known matching verts in the
	same way.

	Example:
	Starting with two meshes and a "matched" vertex selection on two meshes
	I'm calling M1 and M2

	Say:
	v6 on M1 is an edge away from (v3 ,v5) and a face away from (v3 ,v5,v4)
	v9 on M2 is an edge away from (v13,v2) and a face away from (v13,v2,v6)
	And we know going in that: M1.v3=M2.v13, M1.v5=M2.v2, M1.v4=M2.v6

	Then we can say M1.v6=M2.v9 becaue if we substitute all of our known
	matches, we can see that the two vertices are equivalent

	Args:
		M1: A Mesh object
		M2: A second Mesh object
		vertPairs: A List of 2-Tuples of vertex indices
			that are known matches
		matchedNum: The total number of vertices matched up to this point
		vertNum: The total number of vertices in the mesh, for
			percentageDone purposes. Default: None
		symmetry: Boolean value indicating whether the vertPairs
			are mirrored indices on the same mesh. Default: False

	Returns:
		A list of (Vertex,Vertex) pairs that define 1:1 matches

	"""
	orderVerts, shapeVerts = zip(*vertexPairs)
	orderVertsGrow = set(orderVerts)
	shapeVertsGrow = set(shapeVerts)
	orderVerts = set(orderVerts)
	shapeVerts = set(shapeVerts)
	centerVerts = set()

	orderToShape = {s:t for s, t in vertexPairs}

	if vertNum is not None:
		vertNum = float(vertNum) #just for percentage reporting
		if symmetry:
			vertNum /= 2

	if not matchedNum:
		matchedNum = 0

	updated = True
	#counter = 0
	while updated:
		updated = False
		if vertNum is not None:
			percent = (len(orderVerts) + matchedNum) / vertNum * 100
			if pBar is None:
				print "\rPercentage processed: {0:.2f}%".format(percent),
			else:
				pBar.setValue(int(percent))
				QApplication.processEvents()

		# Grow the vert selection along edges and save as a set
		# Grow the vert selection along faces and save as another set
		# Remove already selected verts from both lists
		#
		# The dicts are structured like this:
		#	 {vertex: (orderVertTuple), ...}

		if symmetry:
			# A fun little hack that allows me to treat the left and
			# right hand sides of a model with symmetrical topology as
			# two different meshes to match
			allSet = orderVerts|shapeVerts|centerVerts
			lAllSet = allSet
			rAllSet = allSet
		else:
			lAllSet = orderVerts
			rAllSet = shapeVerts

		orderEdgeDict, orderFaceDict, orderVertsGrow = growTracked(orderMesh, orderVertsGrow, lAllSet)
		shapeEdgeDict, shapeFaceDict, shapeVertsGrow = growTracked(shapeMesh, shapeVertsGrow, rAllSet)

		# if a key has a *unique* (face & edge) value
		# we can match it to the shape
		#
		# so flip the dicts and key off of *BOTH* values
		# simultaneously
		orderEFDict = flipMultiDict(orderEdgeDict, orderFaceDict)
		shapeEFDict = flipMultiDict(shapeEdgeDict, shapeFaceDict)

		# Then, if the swapped dict's value only has 1 item
		# it is a uniquely identified vertex and can be matched
		for orderKey, orderValue in orderEFDict.iteritems():
			if len(orderValue) == 1:
				edgeKey = frozenset(orderToShape[i] for i in orderKey[0])
				faceKey = frozenset(orderToShape[i] for i in orderKey[1])

				try:
					shapeValue = shapeEFDict[(edgeKey, faceKey)]
				except KeyError:
					area = list(edgeKey) + list(faceKey)
					area = list(set(area))
					vp = vertexPairs[:]
					iidx, jidx = zip(*vp)
					m = 'Order {0} to Shape {1}: Match produced no results: Check this order area {2}'.format(iidx, jidx, area)
					#if vertNum is not None and pBar is None:
						#print #clear the percentage value
					raise TopologyMismatch(m)

				if len(shapeValue) != 1:
					#if vertNum is not None and pBar is None:
						#print #clear the percentage value
					raise TopologyMismatch('Match produced multiple results')

				orderVert = orderValue[0]
				shapeVert = shapeValue[0]
				if (shapeVert == orderVert) and symmetry:
					centerVerts.add(shapeVert)
				else:
					orderToShape[orderVert] = shapeVert
					orderVerts.add(orderVert)
					orderVertsGrow.add(orderVert)
					shapeVerts.add(shapeVert)
					shapeVertsGrow.add(shapeVert)

				#pair = (orderVert, shapeVert)
				updated = True

	if vertNum is not None and pBar is None:
		print #clear the percentage value

	#if vertNum is not None:
		#print "\n"
	return [(k, v) for k, v in orderToShape.iteritems()]


