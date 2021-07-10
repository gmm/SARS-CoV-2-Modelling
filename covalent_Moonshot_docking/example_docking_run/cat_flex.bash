#!/bin/bash

for f in *_flex_noheadnotail.pdbqt
do cat flex_header.pdbqt $f flex_tail.pdbqt > ''${f::-22}'flex.pdbqt'
done	
