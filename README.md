# bst-duneuro
**Integration between [Brainstorm](https://neuroimage.usc.edu/brainstorm/) and [Duneuro](http://www.duneuro.org/) projects.**

This project is a spin-off from brainstorm3, Brainstorm Toolbox's main repository.
It is intended to be a place to store all the necesary tools and configuration files in order to recompile **bst_duneuro**, the application which allows Brainstorm to compute FEM-based head models.

###### You can read about this whole integration here 
https://github.com/brainstorm-tools/brainstorm3/issues/242

The first and full discussion regarding the FEM integration :

https://github.com/brainstorm-tools/brainstorm3/issues/185

## Download
Distribution within Brainstorm:
- This is available for download on the Brainstorm website (https://neuroimage.usc.edu/bst/download.php / bst_duneuro_YYMMDD.zip)
- The .zip package is created by calling manually the script: http://neuroimage.usc.edu/bst/update_duneuro.php?g=p845Jui
- It can be downloaded directly with the URL: https://neuroimage.usc.edu/bst/getupdate.php?d=bst_duneuro_200110.zip
- To get the current online version of the package (to know whether it should be updated in Brainstorm): https://neuroimage.usc.edu/bst/getversion_duneuro.php

## Development

Duneuro is based on Dune project, developed in C++ for linux based OS. Thus, there is no windows version directly available. Windows binary application is built by cross-compiling with cygwin. In order to re-compile the application you will need a linux OS and follow these steps:

1. Update needed tools
   in a new 18.04 lts ubuntu we need the following packages:
   ```
   sudo apt-get install git pkg-config cmake mingw-w64 g++-mingw-w64 libc6-dev-i386 python3-pip libeigen3-dev
   ```
   
2. Update the system
  ```
  sudo apt-get install update
  ```
3. clone this repository
```
git clone [THIS REPO'S URL]
```
