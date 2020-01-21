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
 In a new 18.04 lts ubuntu we need the following packages:
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
4. navigate to the repository's directory ```cd bst-duneuro``` and execute compilation tool 

```
chmod +x config/setup_bst_duneuro.sh
config/setup_bst_duneuro.sh
```

The setup script ```setup_bst_duneuro.sh``` is designed to help through the compilation process even if the user's experience with C++, linux and development in general is low. Basically you can execute ```config/setup_bst_duneuro.sh all``` and it will download, configure, compile and serve the binary files. This whole chain has proven to work on a freshly downloaded Ubuntu instance.

In case further development has to be made, you can ```clean``` previous builds and ```download``` the source code before compiling anything. You can also compile only one operating system by calling ```config/setup_bst_duneuro.sh build windows``` or ```linux```. And if you want to test small modifications to your application and you don't want to compile the whole duneuro project in every iteration you can use ```config/setup_bst_duneuro.sh rebuild windows duneuro-matlab```.

See ```config/setup_bst_duneuro.sh``` with no input to read see more examples.

