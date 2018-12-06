#pylint:disable=missing-docstring,invalid-name
'''
Here are a bunch of functions, generators, and co-processes
that are made to allow for a UI easily increment over a set
of possibilities defined by my heuristics.

'''
import itertools
from MeshCrawler.meshcrawlerLib import matchByTopology
from MeshCrawler.meshcrawlerErrors import TopologyMismatch, IslandMismatch

from Qt.QtWidgets import QApplication

###################################################
###				  Approximation 				###
###################################################

def unscrambleMeshByDistance(clean, dirty):
	cleanVerts = clean.vertArray
	dirtyVerts = dirty.vertArray
	return unscrambleByDistance(cleanVerts, dirtyVerts)

def unscrambleByDistance(cleanVerts, dirtyVerts):
	'''
	Given two meshes whose vertices are generally close to one another,
	find a 1:1 mapping where the distances between the mappings are
	minimized.
	This uses the Munkres (aka Hungarian) algorithm and it will *not* map
	more than one vertex to any other

	This is an O(n**3) algorithm, so this is gonna be SLOW for big meshes
	'''
	from scipy.optimize import linear_sum_assignment
	from scipy.spatial.distance import cdist

	dist = cdist(cleanVerts, dirtyVerts)
	idxs = linear_sum_assignment(dist)
	# The clean index will be sorted
	return sorted(zip(*idxs))

def unscrambleByDistance_Pure(cleanVerts, dirtyVerts, invert=False):
	# use the pure python implementation
	# because scipy isn't easy to get for Maya :(
	from munkres import Munkres
	m = Munkres()
	costs = buildCosts(cleanVerts, dirtyVerts, invert=invert)
	indexes = m.compute(costs)
	return sorted(indexes)

def buildCosts(orderCenters, shapeCenters, invert=False):
	# Builds a preference list based on the distance
	# between bounding box centers
	squaredDistances = []
	mul = -1 if invert else 1
	for oC in orderCenters:
		row = []
		for sC in shapeCenters:
			row.append(mul * sum((i-j)**2 for i, j in zip(oC, sC)))
		squaredDistances.append(row)
	return squaredDistances


###################################################
###			 Find pairs from positions			###
###################################################

def _getMinListSizeKey(d):
	lenDict = {}
	for k, v in d.iteritems():
		lenDict.setdefault(len(v), set()).add(k)
	minLen = min(lenDict.iterkeys())
	return lenDict[minLen]

def _getValenceDict(mesh, verts):
	meshValence = {}
	for vert in verts:
		valence = len(mesh.vertNeighbors[vert])
		meshValence.setdefault(valence, []).append(vert)

	return meshValence

def _getMinValencePoints(order, shape, orderVerts, shapeVerts):
	# get the smallest group of common valence points
	orderValence = _getValenceDict(order, orderVerts)
	shapeValence = _getValenceDict(shape, shapeVerts)

	if len(orderValence) != len(shapeValence):
		ovk = set(orderValence.keys())
		svk = set(shapeValence.keys())
		oCheck = []
		sCheck = []
		for key in ovk - svk:
			oCheck.extend(orderValence[key])
		for key in svk - ovk:
			sCheck.extend(shapeValence[key])

		print "Raising Topo Mismatch", oCheck, sCheck
		raise TopologyMismatch("Valence Points Mismatch. Check Order Here {0}, and Shape Here {1}".format(oCheck, sCheck))
	else:
		for key in orderValence:
			if len(orderValence[key]) != len(shapeValence[key]):
				if len(orderValence[key]) < 10 and len(shapeValence[key]) < 10:
					oCheck = [i for i in orderValence[key]]
					sCheck = [i for i in shapeValence[key]]
					raise TopologyMismatch("Valence Points Mismatch. Check Order Here {0}, and Shape Here {1}".format(oCheck, sCheck))
				else:
					raise TopologyMismatch("Valence Points Mismatch. Too many mismatches to be useful")

	minValence = _getMinListSizeKey(orderValence)
	minValence = minValence.pop()

	orderPoints = orderValence[minValence]
	shapePoints = shapeValence[minValence]

	return minValence, orderPoints, shapePoints

def _getNearestGrow(mesh, points, valence):
	# make sure only to run this if there's more than 1 point
	growLength = {} # steps:[vertList]
	for point in points:
		grown = set([point])
		exclude = set()
		steps = 0
		found = False
		while not found and len(grown) > 0:
			grown, exclude = _growByEdge(mesh, grown, exclude)
			steps += 1
			for g in grown:
				if len(mesh.adjacentVertsByEdge(g)) == valence:
					growLength.setdefault(steps, []).append(point)
					found = True
					break
		if not found:
			raise TopologyMismatch("Could not find any other valence points")
	return growLength

def _growByEdge(mesh, growSet, exclude):
	""" Grow a set of verts along edges without any fanciness
	Args:
		mesh: The mesh object for the growth
		growSet: A set of Vertex objects to grow.
		exclude: A set of Vertex objects to exclude from
			the growth

	Returns:
		newGrowSet: the grown verts
		newExclude: the next set of verts to exclude
	"""
	grown = set()
	for vert in growSet:
		grown.update(mesh.adjacentVertsByEdge(vert))

	newgrown = grown - exclude
	newexclude = exclude | growSet

	return newgrown, newexclude

def _growByFace(mesh, growSet, exclude):
	""" Grow a set of verts along edges without any fanciness
	Args:
		mesh: The mesh object for the growth
		growSet: A set of Vertex objects to grow.
		exclude: A set of Vertex objects to exclude from
			the growth

	Returns:
		newGrowSet: the grown verts
		newExclude: the next set of verts to exclude
	"""

	grown = set()
	for vert in growSet:
		grown.update(mesh.adjacentVertsByFace(vert))

	newgrown = grown - exclude
	newexclude = exclude | growSet

	return newgrown, newexclude

def findPossiblePairsByValenceSteps(order, shape, orderVerts, shapeVerts):
	"""
	Makes a list of the vertices that have a specific valence
		(specifically, the valence with the lowest number of verts)
	Then matches the vertices by finding the "grow distance" to the
	closest vertex of the same valence.
	Example:
		There are 54 valence 3 vertices
		There are 13112 valence 4 vertices
		There are 28 valence 5 vertices
		So we use valence 5 vertices

		On the orderMesh, loop through the valence 5 vertices
		Pt 324 is 6 growIterations away from another valence 5 vertex
		Pt 10545 is 6 growIterations away from another valence 5 vertex
		Pt 1484 is 5 growIterations away from another valence 5 vertex
		... and so on

		On the shapeMesh, loop through the valence 5 vertices
		Pt 575 is 6 growIterations away from another valence 5 vertex
		Pt 12245 is 6 growIterations away from another valence 5 vertex
		Pt 3177 is 5 growIterations away from another valence 5 vertex
		... and so on

		There is only one that has a minimum of 5 grows to another
		valence 5 vertex, therefore
		orderMesh.vertices[1484] pairs with shapeMesh.vertices[3177]
		I should be OK if I find a group of 5 or fewer
	"""
	try:
		valence, orderPoints, shapePoints = _getMinValencePoints(order, shape, orderVerts, shapeVerts)
	except KeyError:
		return [], []
	
	if len(orderPoints) == 1:
		return orderPoints, shapePoints

	orderSteps = _getNearestGrow(order, orderPoints, valence)
	shapeSteps = _getNearestGrow(shape, shapePoints, valence)

	orderMinKey = _getMinListSizeKey(orderSteps)
	shapeMinKey = _getMinListSizeKey(shapeSteps)

	common = orderMinKey & shapeMinKey
	if not common:
		return [], []
	minKey = common.pop()

	orderMatches = orderSteps[minKey]
	shapeMatches = shapeSteps[minKey]

	return orderMatches, shapeMatches

def partitionIslands(mesh):
	allverts = set(xrange(len(mesh.vertArray)))
	islands = []
	while allverts:
		seed = set([allverts.pop()])
		island = set()
		while seed:
			seed, island = _growByEdge(mesh, seed, island)
		islands.append(island)
		allverts.difference_update(island)
	return islands

def bbCenter(mesh, island):
	verts = [mesh.vertArray[i] for i in island]
	xAxis, yAxis, zAxis = zip(*verts)
	lC = (min(xAxis), min(yAxis), min(zAxis)) # lowerCorner
	uC = (max(xAxis), max(yAxis), max(zAxis)) # upperCorner
	center = [(i+j)/2.0 for i, j in zip(lC, uC)]
	return center

def makeIslandMarriages(orderMesh, shapeMesh, orderIslands, shapeIslands):
	if len(orderIslands) == 1:
		return [0]

	orderCenters = [bbCenter(orderMesh, i) for i in orderIslands]
	shapeCenters = [bbCenter(shapeMesh, i) for i in shapeIslands]
	unscrambled = unscrambleByDistance_Pure(orderCenters, shapeCenters)
	return [i[1] for i in unscrambled]

def getIslandFaceCount(mesh, island):
	faceSet = set()
	for v in island:
		faceSet.update(mesh.vertToFaces[v])
	return len(faceSet)


###################################################
###              Match CoProcesses              ###
###################################################

def axisMatchGenerator(order, shape, orderIsland, shapeIsland, deep=False):
	''' Find equal valence points that could be matches '''
	orderMatches, shapeMatches = findPossiblePairsByValenceSteps(order, shape, orderIsland, shapeIsland)
	if not orderMatches:
		return
	# There should be an equal number of order and shape matches
	# Build a distance-weighted pairing to minimize
	# the chance of flipping symmetrical meshes

	if len(shapeMatches) < 30:
		orderPoints = [order.vertArray[i] for i in orderMatches]
		shapePoints = [shape.vertArray[i] for i in shapeMatches]
		pairs = unscrambleByDistance_Pure(orderPoints, shapePoints)
		orderIdxs, shapeIdxs = zip(*pairs)
		orderMatches = [orderMatches[i] for i in orderIdxs]
		shapeMatches = [shapeMatches[i] for i in shapeIdxs]

	if not deep:
		for orderPoint in orderMatches:
			for shapePoint in shapeMatches:
				yield (orderPoint, shapePoint)
	else:
		for perm in itertools.permutations(shapeMatches):
			yield list(zip(orderMatches, perm))

def starMatchGenerator(order, shape, orderPoint, shapePoint, reverse=False):
	''' Set up 3-vert neighbors around an axis for an actual match attempt '''
	orderVerts = order.adjacentVertsByEdge(orderPoint)[:2] + [orderPoint]
	shapeStar = shape.adjacentVertsByEdge(shapePoint)
	for i in range(len(shapeStar)):
		#shapeVerts = shapeStar[i:] + shapeStar[:i] + [shapePoint]
		rotStar = shapeStar[i:] + shapeStar[:i]
		if not reverse:
			shapeVerts = rotStar[:2] + [shapePoint]
		else:
			shapeVerts = rotStar[:2][::-1] + [shapePoint]

		yield zip(orderVerts, shapeVerts)

def fullMatchGenerator(order, shape, orderIsland, shapeIsland):
	''' Set up 3-vert neighbors around an axis for an actual match attempt '''
	for orderPoint, shapePoint in axisMatchGenerator(order, shape, orderIsland, shapeIsland):
		for match in starMatchGenerator(order, shape, orderPoint, shapePoint):
			yield match

def fullMatchCoProcess(orderMesh, shapeMesh, skipMismatchedIslands=False):
	''' A co process to iterate most efficiently over all the possible matches for auto-crawling '''
	img = islandMatchCoProcess(orderMesh, shapeMesh, skipMismatchedIslands)
	try:
		oi, si, islandNum = next(img)
		while True:
			smg = fullMatchGenerator(orderMesh, shapeMesh, oi, si)
			found = False
			for sm in smg:
				found = (yield sm, islandNum)
				if found:
					break
			oi, si, islandNum = img.send(found)
	except StopIteration:
		pass

def autoCrawlMeshes(orderMesh, shapeMesh, skipMismatchedIslands=False, pBar=None):
	"""
	Crawl both the order and shape meshes using my heuristics to find
	any matching islands
	"""
	if pBar is not None:
		pBar.setLabelText("Finding Islands")
		pBar.setValue(66)
		QApplication.processEvents()

	fmg = fullMatchCoProcess(orderMesh, shapeMesh, skipMismatchedIslands=skipMismatchedIslands)
	matches = []
	errors = {}
	matchCount = 0
	check = 0
	found = False
	try:
		# Execute the generator up to the first yield, and get the data from it
		sm, curIdx = fmg.send(None) #Send nothing the first time
		idxErrors = []
		while True:
			if pBar is not None:
				pBar.setLabelText("Crawling iteration {0}".format(check))
				pBar.setValue(0)
				QApplication.processEvents()
			check += 1
			try:
				print
				print "Checking Vertex Match", zip(*sm)
				match = matchByTopology(orderMesh, shapeMesh, sm,
					matchedNum=matchCount, vertNum=len(orderMesh.vertArray), pBar=pBar)
			except TopologyMismatch as err:
				idxErrors.append(str(err))
			else:
				matches.append(match)
				found = True
				matchCount += len(match)
			# send the return value into the generator,
			# execute up until the next yield and get the data from it
			sm, idx = fmg.send(found)
			if curIdx != idx:
				if not found:
					errors[curIdx] = idxErrors
				idxErrors = []
				curIdx = idx
			found = False
	except StopIteration:
		if not found:
			raise TopologyMismatch("No Match Found")
	return matches

def islandMatchCoProcess(orderMesh, shapeMesh, deep=False, skipMismatchedIslands=False):
	''' Build possible island matches.
	Either:
	All at once with no user intervention (deep=True) or
	One at a time, waiting on whether it produces a match (deep=False)
	'''
	orderIslandSets = partitionIslands(orderMesh)
	shapeIslandSets = partitionIslands(shapeMesh)

	oiDict = {}
	for oi in orderIslandSets:
		pointCount = len(oi)
		faceCount = getIslandFaceCount(orderMesh, oi)
		oiDict.setdefault((pointCount, faceCount), []).append(oi)

	siDict = {}
	for si in shapeIslandSets:
		pointCount = len(si)
		faceCount = getIslandFaceCount(shapeMesh, si)
		siDict.setdefault((pointCount, faceCount), []).append(si)

	if (len(siDict) == 1 and len(oiDict) == 1 and
			len(siDict.values()[0]) == 1 and len(oiDict.values()[0]) == 1):
		# shortcut the single-island case
		oi = oiDict.values()[0][0]
		si = siDict.values()[0][0]
		if deep:
			yield [[oi]], [[si]]
		else:
			yield oi, si, 0 # don't care about sent values in this case
	else:
		allKeys = set(oiDict.keys() + siDict.keys())
		badKeys = set([key for key in allKeys if key not in oiDict or key not in siDict])
		if badKeys and not skipMismatchedIslands:
			raise IslandMismatch('There are missing islands with these vert/face counts: {0}'.format(list(badKeys)))
		allKeys = allKeys - badKeys
		if not deep:
			for key in allKeys:
				oIslands = oiDict[key]
				sIslands = siDict[key]
				shapeOrder = makeIslandMarriages(orderMesh, shapeMesh, oIslands, sIslands)
				used = [False] * len(shapeOrder)
				for oIdx, oi in enumerate(oIslands):
					for sIdx in shapeOrder:
						if used[sIdx]:
							continue
						si = sIslands[sIdx]
						found = (yield oi, si, oIdx)
						if found:
							used[sIdx] = True
							break
		else:
			# this will give a full possible island match
			# We need to set priorities for the possible matches
			# Islands with different counts on the meshes are lowest priority
			# Then the rest of the islands will be ordered so that high number has low priority
			# (Meaning buttons are less important than an underbody)

			prioDict = {}
			for key in allKeys:
				oIslands = oiDict[key]
				sIslands = siDict[key]
				if len(oIslands) != len(sIslands):
					prioDict.setdefault(100000000, []).append(key)
				else:
					prioDict.setdefault(len(oIslands), []).append(key)

			prioKeys = sorted(prioDict.keys(), reverse=True)
			oiList, siList = [], []
			for pk in prioKeys:
				for key in prioDict[pk]:
					ois = oiDict[key]
					sis = siDict[key]
					shapeOrder = makeIslandMarriages(orderMesh, shapeMesh, ois, sis)
					sis = [sis[i] for i in shapeOrder]
					oiList.append(ois)
					siList.append(sis)

			for perm in _deepIslandGen(siList):
				yield oiList, perm

def _deepIslandGen(lists):
	''' Recursive algorithm for iterating over lists of permutations '''
	if not lists:
		yield []
	else:
		myList = lists[0]
		others = lists[1:]
		for perm in itertools.permutations(myList):
			for accum in _deepIslandGen(others):
				yield [perm] + accum

def matchGenerator(orderMesh, shapeMesh, skipMismatchedIslands=False):
	""" Generates full sets of matches (one per island) rather than a single match per island
	This allows for more manual intervention
	"""
	img = islandMatchCoProcess(orderMesh, shapeMesh, skipMismatchedIslands=skipMismatchedIslands, deep=True)
	for oiGroups, siGroups in img:
		amgs = []
		for oiGroup, siGroup in zip(oiGroups, siGroups):
			for oi, si in zip(oiGroup, siGroup):
				amgs.append(list(axisMatchGenerator(orderMesh, shapeMesh, oi, si)))

		ranges = [range(len(i)) for i in amgs]
		for idxs in itertools.product(*ranges):
			yield [amg[i] for amg, i in zip(amgs, idxs)]

def symmetryGenerator(mesh, island):
	orderMatches, shapeMatches = findPossiblePairsByValenceSteps(mesh, mesh, island, island)
	if not orderMatches:
		return

	if len(shapeMatches) < 30:
		orderPoints = [mesh.vertArray[i] for i in orderMatches]
		shapePoints = [mesh.vertArray[i] for i in shapeMatches]
		pairs = unscrambleByDistance_Pure(orderPoints, shapePoints, invert=True)
		orderIdxs, shapeIdxs = zip(*pairs)
		orderMatches = [orderMatches[i] for i in orderIdxs]
		shapeMatches = [shapeMatches[i] for i in shapeIdxs]

	for orderPoint in orderMatches:
		for shapePoint in shapeMatches:
			if orderPoint == shapePoint:
				continue
			for y in starMatchGenerator(mesh, mesh, orderPoint, shapePoint, reverse=True):
				yield y



def test():
	from blur3d.api.classes.mesh import Mesh

	orderPath = r'H:\public\tyler\bagel\Head_NE.obj'
	shapePath = r'H:\public\tyler\bagel\Head_EN_Bad1.obj'
	#outPath = r'H:\public\tyler\bagel\crawl.obj'

	print "Loading Order Mesh"
	orderMesh = Mesh.loadObj(orderPath)
	print "Loading Shape Mesh"
	shapeMesh = Mesh.loadObj(shapePath)

	matches = autoCrawlMeshes(orderMesh, shapeMesh, skipMismatchedIslands=True)
	
	print "DONE", len(matches)












	#orderPossible, shapePossible = findPossiblePairsByValenceSteps(headOrder.vertices, headShape.vertices)
	#print "POSSIBLE PAIRS:", orderPossible, shapePossible

	#minValence, orderPossible, shapePossible = _getMinValencePoints(headOrder, headShape)
	#match = matchPossiblePairs(headOrder, headShape, orderPossible, shapePossible)
	#print "Match Found!!!"
	#match = matchByTopology(headOrder, headShape, vertPairs, len(headOrder.vertices))

	#updateVertPairs(headOrder, match)
	#UpdateMesh()(headOrder, outPath)

if __name__ == "__main__":
	test()


