
# script to install KITCAT into a virtual environment
echo ""
echo -n "Install KITCAT into a virtual enivonment? (y/n): "
read proceed
if [ "$proceed" = "y" ]; then
    basedir=$(pwd)
else
    exit
fi

# install virtual env if not existed
echo ""
if hash virtualenv 2> /dev/null; then
    virtualenv --python=$(which python) kitcat_ve
else
    echo -n "Install virtualenv? (y/n):"
    read install
    if [ "$proceed" = "y" ]; then
        pip install virtualenv
        virtualenv --python=$(which python) kitcat_ve
    else
        exit
    fi
fi

# source virtual environment and install requirements
echo ""
source kitcat_ve/bin/activate
pip install -r requirements.txt


# print out message
echo ""
echo "To activate virtualenv, run:" 
echo "source kitcat_ve/bin/activate"
echo ""