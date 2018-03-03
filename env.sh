###############################################################
# Set up ROOT
###############################################################
if [ -z $ROOTSYS ]; then
    echo "ROOTSYS is not set. Check your ROOT installation."
fi

###############################################################
# Choose the directory to hold root files
###############################################################
OUTPUT="root"
for name in $OUTPUT
do
    if [ ! -e $name ]; then
        echo "This is the first time lauching this setup!"
        read -p "Where shall I put analysis result? [Press enter for ./$name]"
        while [ -n $REPLY ] && [ ! -e $REPLY ]; do
            read -p "You must provide a valid directory!"
        done
        if [ -z $REPLY ]; then
            echo "I will put the result under ./$name!"
            mkdir -p $name
        else
            ln -s $REPLY $name
        fi
    fi
done

###############################################################
# Export enviromental variables
###############################################################
export CDCS8WORKING_DIR=$PWD
export PATH="$PWD/BinaryFiles/bin:$PWD/scripts:$PATH"
