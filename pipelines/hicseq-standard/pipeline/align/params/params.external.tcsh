#!/bin/tcsh

source ./inputs/params/params.tcsh

set aligner = external
set genome = `./code/read-sample-sheet.tcsh $sheet $object genome`
