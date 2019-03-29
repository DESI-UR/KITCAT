#!/bin/bash

python_version=`python -c "import sys; print(\"{}.{}\".format(sys.version_info[0], sys.version_info[1]))"`

if [ $python_version == "3.6" ]
then 
    echo "Updating the version file for python $python_version only"
    git config --global user.email "travis@travis-ci.org"
    git config --global user.name "Travis CI"
    git clone https://${GH_TOKEN}@github.com/DESI-UR/KITCAT.git -b dev kitcat-dev
    cd catgen-dev
    python py/KITCAT/versioning.py
    git commit -m '[skip travis] after the successful build, updating the version number' py/KITCAT/_version.py
    git push
fi
