#! /bin/bash

## Control Job
mol2gjf="1"
rungaus=""

############### Step 1 ###################
## Convert mol2 to gjf 
## Build mol name folder
if [ ! -z $mol2gjf ];then
for dir in *.mol
do
python ../mol_gjf.py $dir
bname=`basename $dir .mol`
mkdir $bname
mv $dir ${bname}.gjf $bname
done

############### Step 2 ###################
## Submit gaussian job on HPCC
elif [ ! -z ${rungaus} ];then
for dir in `ls -d */`
do
bname=`basename $dir`
cd $dir
../g09qsub.sh -j ${bname}.gjf -m 2 -p 4
cd ..
done

fi
