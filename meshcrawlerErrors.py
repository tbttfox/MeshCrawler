class Mismatch(Exception):
	''' Base class for my mismatch exceptions '''
	pass

class IslandMismatch(Mismatch):
	''' Raised if there are different numbers of islands '''
	pass

class TopologyMismatch(Mismatch):
	''' Raised if the topology doesn't match for an island '''
	pass

