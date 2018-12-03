#pylint: disable=unused-argument, too-many-locals
"""
Remove a subdivision level from a mesh's topology, and
Move the vertices so a new subdivision will match the
original as closely as possible
"""
import numpy as np
from itertools import chain, izip_longest
from Qt.QtWidgets import QApplication

def mergeCycles(groups):
	"""
	Take a list of ordered items, and sort them so the
	last item of a list matches the first of the next list
	Then return the groups of lists mashed together
	for instance, with two cycles:
		input:   [(1, 2), (11, 12), (3, 1), (10, 11), (2, 3), (12, 10)]
		reorder: [[(1, 2), (2, 3), (3, 1)], [(10, 11), (11, 12), (12, 13)]]
		output:  [[1, 2, 3], [10, 11, 12, 13]]

	Also, return whether the cycles merged form a single closed group
	"""
	groups = [list(g) for g in groups]
	heads = {g[0]:g for g in groups}
	tails = {g[-1]:g for g in groups}

	headGetter = lambda x: heads.get(x[-1])
	headSetter = lambda x, y: x + y[1:]

	tailGetter = lambda x: tails.get(x[0])
	tailSetter = lambda x, y: y + x[1:]

	searches = ((headGetter, headSetter), (tailGetter, tailSetter))

	out = []
	cycles = []
	while groups:
		g = groups.pop()
		del heads[g[0]]
		del tails[g[-1]]

		for getter, setter in searches:
			while True:
				adder = getter(g)
				if adder is None:
					break
				g = setter(g, adder)
				del heads[adder[0]]
				del tails[adder[-1]]
				adder[:] = []
			groups = [x for x in groups if x]

		cycle = False
		if g[0] == g[-1]:
			g.pop()
			cycle = True
		cycles.append(cycle)
		out.append(g)
	return out, cycles

def grow(neigh, verts, exclude):
	""" Grow the vertex set, also keeping track
	of which vertices we can safely ignore for
	the next iteration
	"""
	grown = set()
	growSet = verts - exclude
	for v in growSet:
		grown.update(neigh[v])
	newGrown = grown - exclude
	newExclude = exclude | growSet
	return newGrown, newExclude

def buildHint(island, neigh, borders):
	""" Find star points that are an even number of grows from an edge """
	borders = borders & island
	if not borders:
		# Well ... we don't have any good way of dealing with this
		# Best thing I can do is search for a point with the least
		# number of similar valences, and return that
		d = {}
		for v in island:
			d.setdefault(len(neigh[v]), []).append(v)

		dd = {}
		for k, v in d.iteritems():
			dd.setdefault(len(v), []).append(k)

		mkey = min(dd.keys())
		return d[dd[mkey][0]][0]

	exclude = set()
	while borders:
		borders, exclude = grow(neigh, borders, exclude)
		borders, exclude = grow(neigh, borders, exclude)
		for b in borders:
			if len(neigh[b]) != 4:
				return b
	return None

def partitionIslands(faces, neigh, pBar=None):
	""" Find all groups of connected verts """
	allVerts = set(chain.from_iterable(faces))
	islands = []
	count = float(len(allVerts))
	while allVerts:
		verts = set([allVerts.pop()])
		exclude = set()
		while verts:
			verts, exclude = grow(neigh, verts, exclude)
		islands.append(exclude)
		allVerts.difference_update(exclude)

		if pBar is not None:
			pBar.set("Detecting Islands", (count - len(allVerts)) / count)

	return islands

def buildUnsubdivideHints(faces, neigh, borders, pBar=None):
	""" Get one vertex per island that was part of the original mesh """
	islands = partitionIslands(faces, neigh)
	hints = []
	for i, isle in enumerate(islands):
		if pBar is not None:
			pBar.set("Gathering Island Hints", i/float(len(islands)))
		hints.append(buildHint(isle, neigh, borders))

	hints = [h for h in hints if h is not None]
	return hints

def getFaceCenterDel(faces, eNeigh, hints, pBar=None):
	"""
	Given a list of hint "keeper" points
	Return a list of points that were created at the
	centers of the original faces during a subdivision
	"""
	vertToFaces = {}
	vc = set()
	for i, face in enumerate(faces):
		for f in face:
			vertToFaces.setdefault(f, []).append(i)
			vc.add(f)

	count = float(len(vc))
	centers = set()
	midpoints = set()
	originals = set(hints)
	queue = set(hints)

	i = 0
	fail = False
	while queue:
		cur = queue.pop()
		if cur in midpoints:
			continue

		if pBar is not None:
			pBar.set("Partitioning Vertices", i/count)
			i += 2 # Add 2 because I *shouldn't* get any midpoints

		midpoints.update(eNeigh[cur])
		t = centers if cur in originals else originals

		for f in vertToFaces[cur]:
			nVerts = faces[f]

			if len(nVerts) != 4:
				fail = True
				continue

			curFaceIndex = nVerts.index(cur)

			half = int(len(nVerts) / 2)
			diag = nVerts[curFaceIndex - half]

			isOrig = diag in originals
			isCtr = diag in centers
			if not isOrig and not isCtr:
				t.add(diag)
				queue.add(diag)
			elif (isCtr and t is originals) or (isOrig and t is centers) or (diag in midpoints):
				fail = True

	return centers, fail

def getBorders(faces):
	"""
	Arguments:
		faces ([[vIdx, ...], ...]): A face representation

	Returns:
		set : A set of vertex indexes along the border of the mesh
	"""
	edgePairs = set()
	for face in faces:
		for f in range(len(face)):
			edgePairs.add((face[f], face[f-1]))
	borders = set()
	for ep in edgePairs:
		if (ep[1], ep[0]) not in edgePairs:
			borders.update(ep)
	return borders

def buildEdgeDict(faces):
	"""
	Arguments:
		faces ([[vIdx, ...], ...]): A face representation

	Returns:
		{vIdx: [vIdx, ...]}: A dictionary keyed from a vert index whose
			values are adjacent edges
	"""
	out = {}
	for face in faces:
		for f in range(len(face)):
			ff = out.setdefault(face[f-1], set())
			ff.add(face[f])
			ff.add(face[f-2])
	return out

def buildNeighborDict(faces):
	"""
	Build a structure to ask for edge and face neighboring vertices
	The returned neighbor list starts with an edge neighbor, and
	proceeds counter clockwise, alternating between edge and face neighbors
	Also, while I'm here, grab the border verts

	Arguments:
		faces ([[vIdx, ...], ...]): A face representation

	Returns:
		{vIdx: [[vIdx, ...], ...]}: A dictionary keyed from a vert index whose
			values are ordered cycles or fans
		set(vIdx): A set of vertices that are on the border
	"""
	fanDict = {}
	for face in faces:
		for i in range(len(face)):
			fanDict.setdefault(face[i], []).append(face[i+1:] + face[:i])

	borders = set()
	out = {}
	for k, v in fanDict.iteritems():
		fans, cycles = mergeCycles(v)
		for f, c in zip(fans, cycles):
			if not c:
				borders.update((f[0], f[-1], k))
		out[k] = fans
	return out, borders

def _fanMatch(fan, uFan, dWings):
	""" Twist a single fan so it matches the uFan if it can """
	uIdx = uFan[0]
	for f, fIdx in enumerate(fan):
		dw = dWings.get(fIdx, [])
		if uIdx in dw:
			return fan[f:] + fan[:f]
	return None

def _align(neigh, uNeigh, dWings):
	""" Twist all the neighs so they match the uNeigh """
	out = []
	for uFan in uNeigh:
		for fan in neigh:
			fm = _fanMatch(fan, uFan, dWings)
			if fm is not None:
				out.append(fm)
				break
	return out

def buildLayeredNeighborDicts(faces, uFaces, dWings):
	"""
	Build and align two neighbor dicts
	for both faces and uFaces which guarantees that the
	neighbors at the same index are analogous (go in the same direction)
	"""
	neighDict, borders = buildNeighborDict(faces)
	uNeighDict, uBorders = buildNeighborDict(uFaces)

	assert borders >= uBorders, "Somehow the unsubdivided borders contain different vIdxs"

	for k, uNeigh in uNeighDict.iteritems():
		neighDict[k] = _align(neighDict[k], uNeigh, dWings)

	return neighDict, uNeighDict, borders

def _findOldPositionBorder(faces, uFaces, verts, uVerts, neighDict, uNeighDict, borders, vIdx, computed):
	"""
	This is the case where vIdx is on the mesh border
	Updates uVerts in-place
	Arguments:
		faces ([[vIdx ...], ...]): The subdivided face structure
		uFaces ([[vIdx ...], ...]): The unsubdivided face structure
		verts (np.array): The subdivided vertex positions
		uVerts (np.array): The unsubdivided vertex positions
		neighDict ({vIdx:[[vIdx, ...], ...]}): Dictionary of neighbor "fans"
		uNeighDict ({vIdx:[[vIdx, ...], ...]}): Dictionary of neighbor "fans"
		vIdx (int): The vertex position to check
		computed (set): A set of vIdxs that have been computed
	"""
	nei = neighDict[vIdx][0]
	nei = [i for i in nei if i in borders]
	assert len(nei) == 2, "Found multi border, {}".format(nei)
	uVerts[vIdx] = 2*verts[vIdx] - ((verts[nei[0]] + verts[nei[1]]) / 2)
	computed.add(vIdx)

def _findOldPositionSimple(faces, uFaces, verts, uVerts, neighDict, uNeighDict, vIdx, computed):
	"""
	This is the simple case where vIdx has valence >= 4
	Updates uVerts in-place
	Arguments
		faces ([[vIdx ...], ...]): The subdivided face structure
		uFaces ([[vIdx ...], ...]): The unsubdivided face structure
		verts (np.array): The subdivided vertex positions
		uVerts (np.array): The unsubdivided vertex positions
		neighDict ({vIdx:[[vIdx, ...], ...]}): Dictionary of neighbor "fans"
		uNeighDict ({vIdx:[[vIdx, ...], ...]}): Dictionary of neighbor "fans"
		vIdx (int): The vertex position to check
		computed (set): A set of vIdxs that have been computed
	"""
	neigh = neighDict[vIdx][0]

	e = [p for i, p in enumerate(neigh) if i%2 == 0]
	f = [p for i, p in enumerate(neigh) if i%2 == 1]

	es = verts[e].sum(axis=0)
	fs = verts[f].sum(axis=0)

	n = len(e)
	term1 = verts[vIdx] * (n / (n-3.0))
	term2 = es * (4 / (n*(n-3.0)))
	term3 = fs * (1 / (n*(n-3.0)))
	vk = term1 - term2 + term3

	uVerts[vIdx] = vk
	computed.add(vIdx)

def _findOldPosition3Valence(faces, uFaces, verts, uVerts, neighDict, uNeighDict, vIdx, computed):
	"""
	This is the complex case where vIdx has valence == 3
	Updates uVerts in-place
	Arguments
		faces ([[vIdx ...], ...]): The subdivided face structure
		uFaces ([[vIdx ...], ...]): The unsubdivided face structure
		verts (np.array): The subdivided vertex positions
		uVerts (np.array): The unsubdivided vertex positions
		neighDict ({vIdx:[[vIdx, ...], ...]}): Dictionary of neighbor "fans"
		uNeighDict ({vIdx:[[vIdx, ...], ...]}): Dictionary of neighbor "fans"
		vIdx (int): The vertex position to check
		computed (set): A set of vIdxs that have been computed

	Returns:
		bool: Whether an update happened

	"""
	neigh = neighDict[vIdx][0]
	uNeigh = uNeighDict[vIdx][0]
	eNeigh = [n for i, n in enumerate(neigh) if i%2 == 0]
	fNeigh = [n for i, n in enumerate(neigh) if i%2 == 1]
	ueNeigh = [n for i, n in enumerate(uNeigh) if i%2 == 0]
	#ufNeigh = [n for i, n in enumerate(uNeigh) if i%2 == 1]

	intr = computed.intersection(ueNeigh)
	if intr:
		# Easy valence 3 case. I only need
		# The computed new neighbor
		# The midpoint on the edge to that neighbor
		# The "face" verts neighboring the midpoint

		# Get the matching subbed an unsubbed neighbor indexes
		uNIdx = intr.pop()
		nIdx = eNeigh[ueNeigh.index(uNIdx)]

		# Get the "face" verts next to the subbed neighbor
		xx = neigh.index(nIdx)
		fnIdxs = (neigh[xx-1], neigh[(xx+1)%len(neigh)])

		# Then compute
		#vk = 4*k1e - ke - k1fNs[0] - k1fNs[1]
		vka = uVerts[uNIdx] + verts[fnIdxs[0]] + verts[fnIdxs[1]]
		vkb = verts[nIdx] * 4
		uVerts[vIdx] = vkb - vka
		computed.add(vIdx)
		return True
	else:
		# The Hard valence 3 case. Made even harder
		# because the paper has a mistake in it

		# vk = 4*ejk1 + 4*ejpk1 - fjnk1 - fjpk1 - 6*fjk1 + sum(fik)
		# where k1 means subdivided mesh
		# where j means an index, jn and jp are next/prev adjacents
		# sum(fik) is the sum of all the points of the face that
		#     *aren't* the original, or edge-adjacent
		#     There could be more than 1 if an n-gon was subdivided
		#
		# I wonder: If it was a triangle that was subdivided, what
		# would sum(fik) because there are no verts that fit that
		# description.  I think this is a degenerate case


		# First, find an adjacent face on the unsub mesh that
		# is only missing the neighbors of the vIdx

		fnIdx = None
		fik = None
		fCtrIdx = None
		for x, v in enumerate(fNeigh):
			# working with neigh, but should only ever contain uNeigh indexes
			origFace = set([n for i, n in enumerate(neighDict[v][0]) if i%2 == 1])
			check = (origFace - set(ueNeigh)) - set([vIdx])
			if computed >= check:
				fCtrIdx = v
				fnIdx = x
				fik = sorted(list(check))
				break

		if fnIdx is None:
			# No possiblity found
			return False

		# Then apply the equation from above
		neighIdx = neigh.index(fCtrIdx)
		ejnk1 = neigh[(neighIdx+1)%len(neigh)]
		ejpk1 = neigh[neighIdx-1]
		fjnk1 = fNeigh[(fnIdx+1)%len(fNeigh)]
		fjpk1 = fNeigh[fnIdx-1]
		fjk1 = verts[fCtrIdx]
		sumFik = uVerts[fik].sum(axis=0)
		vk = 4*ejnk1 + 4*ejpk1 - fjnk1 - fjpk1 - 6*fjk1 + sumFik
		uVerts[vIdx] = vk
		computed.add(vIdx)
		return True
	return False

def deleteCenters(meshFaces, uvFaces, centerDel):
	"""
	Delete the given vertices and connected edges from a face representation
	to give a new representation. Also, while we're in the nitty gritty, make
	some connections to keep track of 'ancestors' and 'children'

	Arguments:
		meshFaces ([[vIdx, ...], ...]): Starting mesh representation
		centerDel ([vIdx, ...]): The vert indices to delete.

	Returns:
		[[vIdx, ...], ...]: A new mesh with unsubdivided topology
		{vIdx: [vIdx, ...], ...}: A dict mapping deleted edge vertices
			to their neighboring kept vertices
	"""
	# For each deleted index, grab the neighboring faces,
	# and twist the faces so the deleted index is first
	cds = set(centerDel)
	faceDelDict = {}
	uvDelDict = {}
	uvFaces = uvFaces or []

	for face, uvFace in izip_longest(meshFaces, uvFaces):
		fi = cds.intersection(face)
		# If we are a subdivided mesh, Then each face will have exactly one
		# vertex that is part of the deletion set
		if len(fi) != 1:
			raise ValueError("Found a face with an unrecognized connectivity")
		# Get that one vert
		idx = fi.pop()
		# Each face is a cycle. Rotate the cycle
		# so that idx is first in the list
		rv = face.index(idx)
		rFace = face[rv:] + face[:rv]
		faceDelDict.setdefault(idx, []).append(rFace)

		if uvFace is not None:
			rUVFace = uvFace[rv:] + uvFace[:rv]
			uvDelDict.setdefault(idx, []).append(rUVFace)

	newFaces = []
	nUVFaces = []
	wings = {}
	uvWings = {}
	for idx, rFaces in faceDelDict.iteritems():
		uvFaces = uvDelDict.get(idx, [])
		# The faces are guaranteed to be in a single loop cycle
		# so I don't have to handle any annoying edge cases! Yay!
		faceEnds = {f[1]: (f[2], f[3], uvf) for f, uvf in izip_longest(rFaces, uvFaces)} #face ends

		end = rFaces[-1][-1] # get an arbitrary face to start with
		newFace = []
		nUVFace = []
		while faceEnds:
			try:
				diag, nxt, uvf = faceEnds.pop(end)
			except KeyError:
				print "rFaces", rFaces
				print "fe", faceEnds
				raise
			if uvf is not None:
				nUVFace.append(uvf[2])
				uvWings.setdefault(uvf[1], []).append(uvf[2])
				uvWings.setdefault(uvf[3], []).append(uvf[2])

			newFace.append(diag)
			wings.setdefault(end, []).append(diag)
			wings.setdefault(nxt, []).append(diag)

			end = nxt
		newFaces.append(newFace)
		if nUVFace:
			nUVFaces.append(nUVFace)
	nUVFaces = nUVFaces or None

	return newFaces, nUVFaces, wings, uvWings

def fixVerts(faces, uFaces, verts, neighDict, uNeighDict, borders, pinned, pBar=None):
	"""
	Given the faces, vertex positions, and the point indices that
	were created at the face centers for a subdivision step
	Return the faces and verts of the mesh from before the
	subdivision step. This algorithm doesn't handle UV's yet

	Arguments:
		faces ([[vIdx, ...], ...]): A face topology representation
		verts (np.array): An array of vertex positions
		centerDel ([vIdx, ..]): A list 'face center' vertices

	Returns:
		[[vIdx, ...], ...]: A non-compact face topology representation
		np.array: An array of vertex positions
	"""
	uVerts = verts.copy()
	uIdxs = sorted(list(set([i for i in chain.from_iterable(uFaces)])))

	v3Idxs = []
	# bowtie verts are pinned
	bowTieIdxs = []
	computed = set()
	i = 0
	count = float(len(verts))

	for idx in uIdxs:
		if pBar is not None:
			pBar.set("Solving positions", i/count)
		if len(uNeighDict[idx]) > 1:
			bowTieIdxs.append(idx)
			i += 1
		elif idx in pinned:
			pass
		elif idx in borders:
			_findOldPositionBorder(faces, uFaces, verts, uVerts, neighDict, uNeighDict, borders, idx, computed)
			i += 1
		elif sum(map(len, uNeighDict[idx])) > 6: # if valence > 3
			_findOldPositionSimple(faces, uFaces, verts, uVerts, neighDict, uNeighDict, idx, computed)
			i += 1
		else:
			v3Idxs.append(idx)

	updated = True
	while updated:
		updated = False
		rem = set()
		for idx in v3Idxs:
			up = _findOldPosition3Valence(faces, uFaces, verts, uVerts, neighDict, uNeighDict, idx, computed)
			if not up:
				continue
			if pBar is not None:
				pBar.set("Solving positions", i/count)
			i += 1
			updated = True
			rem.add(idx)
		v3Idxs = list(set(v3Idxs) - rem)

	return uVerts

def getUVPins(faces, borders, uvFaces, uvBorders, pinBorders):
	"""Find which uvBorders are also mesh borders"""
	if uvFaces is None: return set()
	if pinBorders:
		return set(uvBorders)

	pinnit = set()
	for face, uvFace in zip(faces, uvFaces):
		for i in range(len(face)):
			f = face[i]
			pf = face[i-1]

			uv = uvFace[i]
			puv = uvFace[i-1]

			if not(f in borders and pf in borders):
				if uv in uvBorders and puv in uvBorders:
					pinnit.add(puv)
					pinnit.add(uv)

	return uvBorders & pinnit

def collapse(faces, verts, uvFaces, uvs):
	""" Take a mesh representation with unused vertex indices and collapse it """
	vset = sorted(list(set(chain.from_iterable(faces))))
	nVerts = verts[vset]
	vDict = {v: i for i, v in enumerate(vset)}
	nFaces = [[vDict[f] for f in face] for face in faces]

	if uvFaces is not None:
		uvset = sorted(list(set(chain.from_iterable(uvFaces))))
		nUVs = uvs[uvset]
		uvDict = {v: i for i, v in enumerate(uvset)}
		nUVFaces = [[uvDict[f] for f in face] for face in uvFaces]
	else:
		nUVs = None
		nUVFaces = None

	return nFaces, nVerts, nUVFaces, nUVs

def unsubdivide(faces, verts, uvFaces, uvs, hints=None, repositionVerts=True, pinBorders=False):
	"""
	Given a mesh representation (faces and vertices) remove the edges added
	by a subdivision, and optionally reposition the verts

	Arguments:
		faces ([[vIdx ...], ...]): The subdivided face structure
		verts (np.array): The subdivided vertex positions
		hints (None/list/set): An list or set containing vertex indices that
			were part of the original un-subdivided mesh.
			If not provided, it will auto-detect based on topology relative to the border
			If there are no borders, it will pick an arbitrary (but not random) star point
		repositionVerts (bool): Whether or not to calculate the original vert positions
	"""
	eNeigh = buildEdgeDict(faces)
	if hints is None:
		borders = getBorders(faces)
		hints = buildUnsubdivideHints(faces, eNeigh, borders)
	centerDel, fail = getFaceCenterDel(faces, eNeigh, hints)
	assert not fail, "Could not detect subdivided topology with the provided hints"

	uFaces, uUVFaces, dWings, uvDWings = deleteCenters(faces, uvFaces, centerDel)
	if repositionVerts:
		# Handle the verts
		neighDict, uNeighDict, borders = buildLayeredNeighborDicts(faces, uFaces, dWings)
		pinned = set(borders) if pinBorders else []
		uVerts = fixVerts(faces, uFaces, verts, neighDict, uNeighDict, borders, pinned)

		# Handle the UVs
		uvNeighDict, uUVNeighDict, uvBorders = buildLayeredNeighborDicts(uvFaces, uUVFaces, uvDWings)
		uvPinned = getUVPins(faces, borders, uvFaces, uvBorders, pinBorders)
		uUVs = fixVerts(uvFaces, uUVFaces, uvs, uvNeighDict, uUVNeighDict, uvBorders, uvPinned)
	else:
		uVerts = verts
		uUVs = uvs

	rFaces, rVerts, rUVFaces, rUVs = collapse(uFaces, uVerts, uUVFaces, uUVs)
	return rFaces, rVerts, rUVFaces, rUVs

def writeObj(faces, verts, uvFaces, uvs, path):
	""" Write a .obj file """
	import time
	lines = []

	lines.append('# obj export')
	lines.append('# file created: {}'.format(time.strftime('%c')))

	lines.append('')
	lines.append('# begin {} vertices'.format(len(verts)))
	for v in verts:
		lines.append('v {0:.8f} {1:.8f} {2:.8f}'.format(*v))
	lines.append('# end {} vertices'.format(len(verts)))

	if uvs.size > 0:
		lines.append('')
		lines.append('# begin {} uvs'.format(len(uvs)))
		for u in uvs:
			lines.append('vt {0:.8f} {1:.8f}'.format(*u))
		lines.append('# end {} uvs'.format(len(uvs)))

	lines.append('')
	lines.append('# begin {} faces'.format(len(faces)))
	uvFaces = uvFaces or []
	for f, uvf in izip_longest(faces, uvFaces):
		if uvf is not None:
			vnums = ' '.join(['{}/{}'.format(i+1, u+1) for i, u in zip(f, uvf)])
		else:
			vnums = ' '.join([str(i+1) for i in f])
		lines.append('f ' + vnums)

	lines.append('# end {} faces'.format(len(faces)))
	out = '\n'.join(lines)
	with open(path, 'w') as f:
		f.write(out)

def readObj(path):
	""" Read a .obj file """
	vertices = []
	uvs = []
	faces = []
	uvFaces = []
	with open(path, 'r') as inFile:
		lines = inFile.readlines()

	badUvs = False
	for line in lines:
		sp = line.split()
		if sp == []:
			pass
		elif sp[0] == "v":
			vertices.append([float(i) for i in sp[1:4]])
		elif sp[0] == "vt":
			uvs.append([float(i) for i in sp[1:3]])
		elif sp[0] == "f":
			face = []
			for s in sp[1:]:
				vt = [int(i)-1 if i else None for i in s.split('/')]
				# Pad out the face vert/uv/normal triples
				face.append(vt + [None]*3)
			# Then cut them back to 3
			fTrips = [i[:2] for i in face]
			face, uvFace = zip(*fTrips)
			faces.append(face)
			uvFaces.append(uvFace)
			if None in uvFace:
				badUvs = True

	if badUvs:
		uvFaces = None
		uvs = None
	else:
		uvs = np.array(uvs)

	return faces, np.array(vertices), uvFaces, uvs

def test():
	""" Quick test """
	inPath = r'D:\Users\tyler\Desktop\trueUnsub\Subd.obj'
	outPath = r'D:\Users\tyler\Desktop\trueUnsub\UnSub.obj'
	print "Reading"
	faces, verts, uvFaces, uvs = readObj(inPath)
	print "Unsubbing"
	rFaces, rVerts, rUVFaces, rUVs = unsubdivide(faces, verts, uvFaces, uvs)
	print "Writing"
	writeObj(rFaces, rVerts, rUVFaces, rUVs, outPath)

if __name__ == "__main__":
	test()


