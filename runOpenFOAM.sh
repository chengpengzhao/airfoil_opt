#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory
#------------------------------------------------------------------------------
python3 PARSEC.py;

caseName="airfoil_run";
runpath="./cases/${caseName}";
if [ -d ${runpath} ];then
echo "Delete files......\n\n"
rm -rf ${runpath}
fi

foamCloneCase airfoil_template/ ${runpath};
cp blockMeshDict ${runpath}/system;
cd cases/airfoil_run;
blockMesh;
checkMesh;
renumberMesh -overwrite;
paraFoam &
