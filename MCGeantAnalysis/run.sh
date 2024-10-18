#!/bin/bash

#if [ -z ${3+x} ]; then
if [ -z ${1+x} ]; then
   # echo "expected arguments: input_file_name setup_description_json run_number "
   echo "expected arguments: input_file_name "
    exit 
fi

#./MCGeantAnalysis.x -t mcGeant -f $1 -u userParams.json -l $2 -i $3

/home/pooja/framework/build/MCGeantAnalysis/MCGeantAnalysis.x -t mcGeant \
-f /home/pooja/jpet-geant/jpetGeant/build/bin/$1 \
-u /home/pooja/framework/scripts_for_Pooja/userParams.json \
-i 9 \
-l /home/pooja/framework/j-pet-framework-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json \
-o /home/pooja/framework/cat_MCGeantAnalysis/ -b
