#! /bin/bash

for dir in `seq 1 1200`
do
	if [ -d Keyfiles/Delaney_${dir} ]; then
		echo "Doing Delaney_${dir}...."
		python ./Process_Delaney.py Keyfiles/Delaney_${dir}/Delaney_${dir}_resp.pqrta delaney-processed.csv
	else
		echo "Delaney_${dir} not exist.."
	fi
done
