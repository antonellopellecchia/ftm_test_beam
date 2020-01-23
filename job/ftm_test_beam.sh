#!/bin/sh
export HOME="/lustrehome/antonellopellecchia"
source $HOME/sourcefiles/set_geant4.sh
$HOME/simulations/geant4/ftm_test_beam/build/exampleB1 $HOME/simulations/geant4/ftm_test_beam/build/run.mac
