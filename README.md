# bst-duneuro
**Integration between [Brainstorm](https://neuroimage.usc.edu/brainstorm/) and [Duneuro](http://www.duneuro.org/) projects.**

This project is a spin-off from brainstorm3, Brainstorm Toolbox's main repository.
It is intended to be a place to store all the necesary tools and configuration files in order to recompile **bst_duneuro**. This application allows Brainstorm to compute FEM-based head models through duneuro project.

##### You can read about this whole integration here:
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

Duneuro is based on Dune project, developed in C++, in a linux based OS. Thus, there is no windows version directly available. Windows binary application is built by cross-compiling with cygwin. In order to compile the application you will need a linux OS and follow these steps:

1. Update the needed tools/packages and update your system.
In a new 18.04 lts ubuntu you will need the following packages:
```
sudo apt-get install git pkg-config cmake mingw-w64 g++-mingw-w64 libc6-dev-i386 python3-pip libeigen3-dev
sudo apt-get install update
```

2. Clone this repository
```
git clone https://github.com/brainstorm-tools/bst-duneuro.git
```

3. We have developed a setup/configuration semi automatic script called  ```setup_bst_duneuro.sh```. In order to use it, make sure the script is executable and run it.
```
cd bst-duneuro
chmod +x config/setup_bst_duneuro.sh
config/setup_bst_duneuro.sh
```

The setup script is designed to go through the compilation steps even when the user's experience with C++, linux and development in general is scarce. Basically, you can execute ```config/setup_bst_duneuro.sh``` and it will show you different options and examples in its help message. 

This whole chain has proven to work on a fresh/new Ubuntu instance, so we hope you don't find problems with your system, although we are not sure of how previous configurations with MINGW or other sw might possibly interact with this whole build process.

In order to download, configure, compile and copy the binary application files into ```bin/``` execute
```
setup_bst_duneuro.sh all
```

In case further development has to be made, you can ```clean``` previous builds, ```download``` the source code before compiling anything. You can also compile only one operating system by calling ```setup_bst_duneuro.sh build windows``` or ```setup_bst_duneuro.sh build linux```. 

Usual compilations times are about 10-15 minutes for each target operating system. If you want to test small modifications to your application and you don't want to compile the whole duneuro project in every iteration you can use the ```rebuild``` command. This way you will only build the specific part of duneuro needed which interacts directly with bst_duneuro application. This module is ```duneuro-matlab```. You can do:

```
setup_bst_duneuro.sh rebuild windows duneuro-matlab
setup_bst_duneuro.sh rebuild linux duneuro-matlab
```
or
```
setup_bst_duneuro.sh rebuild all duneuro-matlab
```
to rebuld both linux and window's applications

#### Comments on Matlab and mex file implementation. 
If you have Matlab installed in your linux version and you plan to use Brainstorm in this system, you can take advantage of the possibility of building the application as a mex file which will decrease execution times and allow for more interaction between Matlab and duneuro.

In order to do this modify Matlab's path in ```config/config_release_linux.opts``` and modify ```src/duneuro-matlab/toto/CMakeLists.txt``` file in order to add a separate mex application.


