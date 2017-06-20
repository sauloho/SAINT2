#!/bin/bash


# To build side-chains, we recommend using the relax protocol from Rosetta

# Make this variable point to your current Rosetta installation.
export ROSETTA_DIR=/data/icarus/not-backed-up/oliveira/rosetta_src_2016.32.58837_bundle/

$ROSETTA_DIR/main/source/bin/relax.linuxgccrelease -in:file:s $1 -relax:fast
