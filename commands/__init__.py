try:
	from .maxCommands import (cloneObject, createRawObject, freezeObject, getFaces, getObjectByName, getObjectName, getSingleSelection, getUVs, getVerts, getVertSelection, rootWindow, selectVerts, setAllVerts, setObjectName)
except ImportError:
	try:
		from .mayaCommands import (cloneObject, createRawObject, freezeObject, getFaces, getObjectByName, getObjectName, getSingleSelection, getUVs, getVerts, getVertSelection, rootWindow, selectVerts, setAllVerts, setObjectName)
	except ImportError:
		try:
			from .xsiCommands import (cloneObject, createRawObject, freezeObject, getFaces, getObjectByName, getObjectName, getSingleSelection, getUVs, getVerts, getVertSelection, rootWindow, selectVerts, setAllVerts, setObjectName)
		except ImportError:
			from .externalCommands import (cloneObject, createRawObject, freezeObject, getFaces, getObjectByName, getObjectName, getSingleSelection, getUVs, getVerts, getVertSelection, rootWindow, selectVerts, setAllVerts, setObjectName)

