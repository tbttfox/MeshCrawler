try:
	from .maxCommands import (getSingleSelection, getObjectByName, getObjectName,
						  getFaces, getVerts, createRawObject, selectVerts, getVertSelection)
except ImportError:
	try:
		from .mayaCommands import (getSingleSelection, getObjectByName, getObjectName,
						  getFaces, getVerts, createRawObject, selectVerts, getVertSelection)
	except ImportError:
		try:
			from .xsiCommands import (getSingleSelection, getObjectByName, getObjectName,
						  getFaces, getVerts, createRawObject, selectVerts, getVertSelection)
		except ImportError:
			from .externalCommands import (getSingleSelection, getObjectByName, getObjectName,
						  getFaces, getVerts, createRawObject, selectVerts, getVertSelection)

