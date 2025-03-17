#!/bin/bash

#if [ -z ${3+x} ]; then
if [ -z ${1+x} ]; then
   # echo "expected arguments: input_file_name setup_description_json run_number "
   echo "expected arguments: input_file_name "
    exit 
fi

#./MCGeantAnalysis.x -t mcGeant -f $1 -u userParams.json -l $2 -i $3

/data/4/users/ptanty/jpet_geant/framework/build/MCGeantAnalysis/MCGeantAnalysis.x -t mcGeant \
-f /data/4/users/ptanty/jpet_geant/test2/J-PET-geant4/build/bin/$1 \
-u /data/4/users/ptanty/jpet_geant/framework/j-pet-framework-examples/MCGeantAnalysis/userParams.json \
-i 9 \
-l /data/4/users/ptanty/jpet_geant/framework/j-pet-framework-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json \
-o .  -b
