'''
Copyright 2018, Blur Studio

This file is part of MeshCrawler.

MeshCrawler is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MeshCrawler is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with MeshCrawler.  If not, see <http://www.gnu.org/licenses/>.
'''

MESHCRAWLER_UI = None
MESHCRAWLER_UI_ROOT = None
def runMeshCrawlerUI():

	from MeshCrawler.meshcrawlerDialog import MeshCrawlerDialog
	from MeshCrawler.commands import rootWindow

	global MESHCRAWLER_UI
	global MESHCRAWLER_UI_ROOT

	# make and show the UI
	MESHCRAWLER_UI_ROOT = rootWindow()
	# Keep a global reference around, otherwise it gets GC'd
	MESHCRAWLER_UI = MeshCrawlerDialog(parent=MESHCRAWLER_UI_ROOT)
	MESHCRAWLER_UI.show()

