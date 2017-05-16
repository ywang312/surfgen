#!/bin/bash


source $SURFGEN/bin/setsgenvars.sh

echo $SGENFC -I$SURFGEN/source/ getneighbor.f90 -o getneighbor.x $SGENFLAG $SGENLIB
$SGENFC -I$SURFGEN/source/ getneighbor.f90 -o getneighbor.x $SGENFLAG $SGENLIB

echo Done
