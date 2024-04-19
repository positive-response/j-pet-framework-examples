#!/bin/bash

#if [ -z ${3+x} ]; then
if [ -z ${1+x} ]; then
   # echo "expected arguments: input_file_name setup_description_json run_number "
   echo "expected arguments: input_file_name "
    exit 
fi

#./MCGeantAnalysis.x -t mcGeant -f $1 -u userParams.json -l $2 -i $3

#/home/pooja/framework/MCGEANT_build/MCGeantAnalysis/MCGeantAnalysis.x -t mcGeant -f /home/pooja/mc_data/2gamma/$1 -u /home/pooja/framework/j-pet-framework-examples/MCGeantAnalysis/userParams.json -i 9 -l /home/pooja/framework/j-pet-framework-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json -o /home/pooja/framework/data_from_MCAnalysis/2gamma/ -b



/home/pooja/framework/MCGEANT_build/MCGeantAnalysis/MCGeantAnalysis.x -t mcGeant -f $1 -u /home/pooja/framework/j-pet-framework-examples/MCGeantAnalysis/userParams.json -i 9 -l /home/pooja/framework/j-pet-framework-examples/CalibrationFiles/9_RUN/detectorSetupRun9.json -o /home/pooja/framework/data_from_MCAnalysis/5gamma/ -b
