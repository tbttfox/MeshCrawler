# MeshCrawler
Automated Topology Tools

MeshCrawler is a pure python implementation of un-subdivision and topology matching with user facing features that other tools just don't have.

### First and foremost: No Selecting Vertices!!
MeshCrawler has some algorithms to make guesses that will work for 99.99% of cases, including multiple islands for both matching and unsubdivision. But you can select vertices if you need, and the tool will provide you with its guesses and allow you to edit them, so you're not on your own.

### Pure Python (+ numpy)
All of the algorithms are written in pure Python (and numpy) with an abstraction layer. This means as long as you can build a couple python lists and run a Qt ui, you can make MeshCrawler work.

If you're using Maya, you'll have to install a version of numpy to make this work.

## Installation
* Maya Numpy
   * Windows
      * Get numpy for Maya from [here](https://forums.autodesk.com/t5/maya-programming/numpy-1-13-1-scipy-0-19-1-for-maya-2018/td-p/7362541)
      * Change the file name from .whl to .zip
      * Extract the two folders in there (`numpy` and `numpy-1.13.1.dist-info`) into your `Documents\maya\scripts` folder

* Maya MeshCrawler
   * Download the MeshCrawler .zip file from github and extract it. 
   * Rename the folder that contains all the python files from `MeshCrawler-master` to `MeshCrawler`
   * Move that folder into your `Documents\maya\scripts` folder
   * Copy the contents of `shelfBtn.py` into a shelf button

