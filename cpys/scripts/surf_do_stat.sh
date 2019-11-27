#!/bin/bash

path=$1 #'../20_startmods/'
crust=$2 # number of crustal layers, which are supposed to be partly fixed
it=$3 # number of iterations
fix=$4
damp=$5
window=$6
p=$7
ensamble_name=$8


ext='.mod'

for m in $path*.mod
do
	rm tmp*
	core=${m#$path}	
	core=${core%$ext}
	echo "$core"
	cp $m ../start.mod
	INVERT_surf.sh $core $crust $it $fix $damp $window $p 
	#surf_git_stat.sh "$core surf-inverted"
done

nplotter.py -m $ensamble_name *end.mod
#git add models.png
#git commit -m 'models.png'

average.sh $ensamble_name *end.mod
#git add .
#git commit -m 'mean models'
