# MeshCrawler
Automated Topology Tools

MeshCrawler is a pure python implementation of un-subdivision and topology matching with user facing features that other tools just don't have.

### First and foremost: No Selecting Vertices!!
MeshCrawler has some algorithms to make guesses that will work for 99.99% of cases, including multiple islands for both matching and unsubdivision. But you can select vertices if you need, and the tool will provide you with its guesses and allow you to edit them, so you're not on your own.

### Pure Python (+ numpy)
All of the algorithms are written in pure Python (and numpy) with an abstraction layer. This means as long as you can build a couple python lists and run a Qt ui, you can make MeshCrawler work.

Sorry Maya guys, you'll have to install a version of numpy to make this work.

