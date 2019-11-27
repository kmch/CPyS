#!/bin/bash

dir='../20_startmods/' # choice of directory  with ensamble of starting models (path relative to SURF/RF/JOINT!!!)
## obligatory backslash at the end!

# The following must be consistent with INVERT*
window=27 # time window of RF
per=150 # max period for swd [s]
it=20 # number of iterations
nm=20 # number of startmods

#plot.sh -m $dir/*mod

for invtype in SURF RF JOINT
do
cd $invtype
ensamble_name=$invtype$per'_'$nm'nm_'$it'it' # you may add n of start_mods (dir) etc.
echo $ensamble_name

done

########### SURF
cd SURF
ensamble_name='swd'$per'_'$nm'nm_'$it'it' # you may add n of start_mods (dir) etc.
echo $ensamble_name
#
## do statistical inversion and take average
#
#surf_do_stat.sh $ensamble_name $dir
#
## plot mean +/- std model
#
#envelope_plot.sh -em $ensamble_name stdm*.mod stdp*.mod mean*.mod
#
## generate dsp files from all end.mods 
#
# doit_SWD_only.sh *end.mod
#
## plot surf envelope fit
#
#envelope_plot.sh -xd $ensamble_name ../all.dsp *dsp


########### RF
#cd ../RF
#ensamble_name='rf'$per'_'$nm'nm_'$it'it' # you may add n of start_mods (dir) etc.
#echo $ensamble_name
##
### do statistical inversion and take average
##
#rf_do_stat.sh $ensamble_name $dir
##
### plot mean +/- std model
##
#envelope_plot.sh -em $ensamble_name stdm*.mod stdp*.mod mean*.mod
##
### generate a.sac files from all end.mods
##
#doit_RF_only.sh *end.mod
##
### plot rf envelope fit
##
#envelope_plot.sh -er $ensamble_name ../*a.sac *a.sac
##
### show rf fit png
##
##show_rf_fit.sh
##
### plot rf wiggle bestfit
##
#envelope_plot.sh -wr $ensamble_name ../*a.sac 4.00*a.sac # assumes 4.00 gives best fit!!!!!!!!!! (add interactive later)




########### JOINT
cd ../JOINT
ensamble_name='joint'$per'_'$nm'nm_'$it'it' # you may add n of start_mods (dir) etc.
echo $ensamble_name
##
### do statistical inversion and take average
##
#joint_do_stat.sh $ensamble_name $dir
##
### plot mean +/- std model
##
#envelope_plot.sh -em $ensamble_name stdm*.mod stdp*.mod mean*.mod
#
## generate dsp files from all end.mods 
#
doit_SWD_only.sh *end.mod # or doit together with Rf....
#
## plot surf envelope fit
#
envelope_plot.sh -xd $ensamble_name ../all.dsp *dsp
#
## generate a.sac files from all end.mods
#
doit_RF_only.sh *end.mod
#
## plot rf envelope fit
#
envelope_plot.sh -er $ensamble_name ../*a.sac *a.sac
#
## show rf fit png
#
#show_rf_fit.sh
#
## plot rf wiggle bestfit
#
envelope_plot.sh -wr $ensamble_name ../*a.sac 4.00*a.sac # assumes 4.00 gives best fit!!!!!!!!!! (add interactive later)
