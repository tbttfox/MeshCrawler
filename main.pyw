# make sure this is being run as the main process
if ( __name__ in ( '__main__', '__builtin__' ) ):
	
	from MeshCrawler.meshcrawlerDialog import MeshCrawlerDialog
	import blur3d
	blur3d.launch( MeshCrawlerDialog, instance=True )

