#!/bin/bash

doClean () {
  echo -e "\nDeleting previous builds.\n"
  rm -rf build*
  echo -e "\nDeleting previously downloaded code.\n"
  rm -rf src
}

doDownload () {
  echo -e "\nDeleting previously downloaded code.\n"
  rm -rf src
  mkdir src
  cd src

  #core modules
  git clone --branch releases/2.6 https://gitlab.dune-project.org/core/dune-common.git
  git clone --branch releases/2.6 https://gitlab.dune-project.org/core/dune-geometry.git
  git clone --branch releases/2.6 https://gitlab.dune-project.org/core/dune-grid.git
  git clone --branch releases/2.6 https://gitlab.dune-project.org/core/dune-localfunctions.git
  git clone --branch releases/2.6 https://gitlab.dune-project.org/core/dune-istl.git
  #typetree needed
  git clone --branch releases/2.6 https://gitlab.dune-project.org/staging/dune-typetree.git
  git clone --branch releases/2.6 https://gitlab.dune-project.org/staging/dune-uggrid.git
  #extension functions module
  git clone https://gitlab.dune-project.org/staging/dune-functions.git
  cd dune-functions
  git checkout bd847eb9f6617b116f5d6cb4930e5417d7b6e9a7
  git checkout -b c-interface
  cd ..
  #pdelab and  required
  git clone --branch releases/2.6 https://gitlab.dune-project.org/pdelab/dune-pdelab.git
  #duneuro 
  git clone --branch releases/2.6 https://gitlab.dune-project.org/duneuro/duneuro.git 
  git clone --branch feature/c-interface https://gitlab.dune-project.org/duneuro/duneuro-matlab.git
  
  #git clone --branch 2.6-changes https://gitlab.dune-project.org/duneuro/duneuro-tests.git
  #git clone --branch releases/2.6 https://gitlab.dune-project.org/quality/dune-testtools.git
  
  cd ..
}

doConfigure () {

  # Fortran patch
  echo -e "\n\n#####################################################\n\n               Patching fortran issue!!\n\n#####################################################\n"
  sed -i 's/workaround_9220(Fortran Fortran_Works)/if(ENABLE_Fortran)\n    workaround_9220(Fortran Fortran_Works)\n  endif()/g' src/dune-common/cmake/modules/DuneMacros.cmake 

  # Copy C++ files to duneuro-matlab
  echo -e "\nCopying brainstorm_app folder to duneuro-matlab.\n"
  cp -r config/brainstorm_app src/duneuro-matlab

  # Modify CMakeLists file
  if ! grep -Fxq "brainstorm_app" src/duneuro-matlab/CMakeLists.txt ; then
    echo -e "\nAdding subdirectory brainstorm_app to cmake lists file.\n"
    sed -i 's#add_subdirectory("cmake/modules")#add_subdirectory("cmake/modules")\nadd_subdirectory("brainstorm_app")#g' src/duneuro-matlab/CMakeLists.txt
  fi

  # Modify full path to toolchain file in windows opts file.
  echo -e "\nAdding full path to toolchain file.\n"
  sed "s~<<fullpath>>~${BASE_DIR}~g" config/config_release_windows.opts.template > config/config_release_windows.opts
}

doBuild () {
  echo -e "\nDeleting previous builds.\n"
  rm -rf ${BUILD_DIR}*

  echo -e "\nBuilding duneuro. Log file created in ${BUILD_DIR}_log.txt\n"
  # echo `pwd`
  # echo config file = $CONFIG_FILE
  # echo build dir = $BUILD_DIR
  src/dune-common/bin/dunecontrol --opts=config/${CONFIG_FILE} --builddir=`pwd`/${BUILD_DIR} all | tee "${BUILD_DIR}_log.txt"
 
}

doReBuild () {
  echo -e "\nBuilding duneuro. Log file created in ${BUILD_DIR}_rebuild_log.txt\n"
  src/dune-common/bin/dunecontrol --opts=config/${CONFIG_FILE} --builddir=`pwd`/${BUILD_DIR} --only=${MODU} all | tee "${BUILD_DIR}_rebuild_log.txt"
}

doHandleFiles () {
  
  if test -f "${BUILD_DIR}/duneuro-matlab/brainstorm_app/bst_duneuro${EXTENSION}"; then
    echo -e "\n Everything went fine!\n"
    mv ${BUILD_DIR}/duneuro-matlab/brainstorm_app/bst_duneuro${EXTENSION} bin/bst_duneuro_$(date +"%d_%m_%Y")${EXTENSION}
    echo -e "\n Brainstorm - Duneuro application should be in /bin\n"
    else
    echo -e "\n Something went wrong. Check ${BUILD_DIR}_log.txt or ${BUILD_DIR}_rebuild_log.txt .\n"
  fi
}

print_style () {

    if [ "$2" == "info" ] ; then
        COLOR="96m";
    elif [ "$2" == "success" ] ; then
        COLOR="92m";
    elif [ "$2" == "warning" ] ; then
        COLOR="93m";
    elif [ "$2" == "danger" ] ; then
        COLOR="91m";
    else #default color
        COLOR="0m";
    fi

    STARTCOLOR="\e[$COLOR";
    ENDCOLOR="\e[0m";

    printf "$STARTCOLOR%b$ENDCOLOR" "$1";
}

doPrintHelp () {
  print_style "\n  Duneuro compilation helper\n" "warning";
  print_style "\n" ;
  print_style "  If everything goes fine, an application will be copied to wherever\n" ;
  print_style "  this script is called from.\n" ;
  print_style "\n" ;
  print_style "  Usage: setup_bst_duneuro.sh [task] [vars...]\n" "warning";
  print_style "\n" ;
  print_style "  Tasks: \n" ;
  print_style "      - clean    clean all downloaded code and previous builds\n" ;
  print_style "      - download download source code\n" ;
  print_style "      - build  [OS]  build duneuro\n" ;
  print_style "                 This will also delete previous builds and create\n" ;
  print_style "                 a log txt file. Default = windows\n" ;
  print_style "      - rebuild [OS] [MODULE]  rebuild only one module\n" ;
  print_style "                 This will also create a log txt file.\n" ;
  print_style "                 Default OS = windows . Default MODULE = duneuro-matlab \n" ;
  print_style "  Operating System Versions: \n" ;
  print_style "      - windows  windows version.\n" ;
  print_style "      - linux    linux and mac version.\n" ;
  print_style "\n" ;
  print_style "  Examples: \n" ;
  print_style "      - setup_bst_duneuro all            Delete previous builds and downloads.\n" "warning";
  print_style "                                           Then build all versions\n" "warning";
  print_style "      - setup_bst_duneuro.sh clean       Delete all previous builds and source\n" ;
  print_style "      - setup_bst_duneuro.sh download    Download the necesary repositories.\n" ;
  print_style "      - setup_bst_duneuro build windows  Build windows version\n" ;
  print_style "      - setup_bst_duneuro build linux    Build linux version\n" ;
  print_style "      - setup_bst_duneuro rebuild windows matlab-duneuro\n" ;
  print_style "                                         Rebuild only matlab-duneuro module in its\n";
  print_style "                                           windows version.\n";
  print_style "      - setup_bst_duneuro rebuild linux matlab-duneuro  \n" ;
  print_style "                                         Rebuild only matlab-duneuro in its windows\n" ;
  print_style "                                           version.\n" ;
  print_style "\n\n\n" ;
}

setVariables () {
  BUILD_DIR="build_release_${OS}"
  CONFIG_FILE=config_release_${OS}.opts
  if [[ $OS == windows ]]; then
    EXTENSION=".exe"
  elif [[ $OS == linux ]]; then
    EXTENSION=""
  fi
}

#start script
CURRENT_DIR=`pwd`
cd $(dirname "$0")/..
BASE_DIR=`pwd`
SCRIPT_NAME=$(basename "$0")
#echo Running from ${BASE_DIR} ...


if [[ -z "$1" ]]; then
  doPrintHelp
  exit 1
elif [[ $1 == clean ]]; then
  doClean
elif [[ $1 == download ]]; then
  doDownload
  doConfigure
elif [[ $1 == build ]]; then
  if [[ -z "$2" ]]; then
    doPrintHelp
    exit 1
  elif [[ $2 == windows || $2 == linux ]]; then
    OS=$2
    setVariables
    doBuild
    doHandleFiles
  elif [[ $2 == all ]]; then
    config/${SCRIPT_NAME} build windows
    config/${SCRIPT_NAME} build linux
  else
    doPrintHelp
    exit 1
  fi
elif [[ $1 = all ]]; then
  config/${SCRIPT_NAME} clean
  config/${SCRIPT_NAME} download
  config/${SCRIPT_NAME} build all
elif [[ $1 = rebuild ]]; then
  if [[ -z "$2" ]]; then
    doPrintHelp
    exit 1
  elif [[ $2 == windows || $2 == linux ]]; then
    OS=$2
    if [[ -z "$3" ]]; then
      doPrintHelp
      exit 1
    else
      MODU=$3
      setVariables
      doReBuild
      doHandleFiles
    fi
  elif [[ $2 == all ]]; then
    config/${SCRIPT_NAME} rebuild windows $3
    config/${SCRIPT_NAME} rebuild linux $3
  else
    doPrintHelp
    exit 1
  fi
else
  doPrintHelp
  exit 1
fi

exit 0 #SUCCESS

