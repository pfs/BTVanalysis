#!/bin/bash

outDir="/afs/cern.ch/user/p/psilva/public/html/BTV/7TeV_legacy/"
mkdir -p ${outDir}
wps=("0" "2")
for i in ${wps[@]}; do

    #prepare templates and build datacards
    #root -b -q "drawObservedDiffBtagDistribution.C+(\"legacy7TeV_bdt/plotter.root\",\"csv${i}\")";
    #root -b -q "buildDataCardsForDifferentialMeasurement.C+(\"csv${i}/csv${i}_btags.root\",\"legacy7TeV_bdt/plotter.root\",\"legacy7TeV_bdt/plotter_syst.root\",\"csv${i}\")";
    #sed -e 's/_TAGGER_/csv/g' -e 's/_WP_/'${i}'/g' < diffbtag_index.html.template > #{outDir}/ftm-results/csv${i}/index.html

    #move to output directory
    #mv -ir diffbtag ${outDir}/ftm-results
    
    #run fit
    rundFtM --flav ${outDir}/ftm-results/csv${i}/flavbreakup.json \
	--btag ${outDir}/ftm-results/csv${i}/csv${i}_eff.json \
	--in ${outDir}/ftm-results/csv${i}/csv${i}_btags.root --fit 0 --fix 1,5,6 --inc;
done

#a=(`ls ${outDir}/ftm-results/csv*/*workspace.root`)
#allWS=""
#for i in ${a[@]}; do
#    if [ -z ${allWS} ]; then
#	allWS=${i};
#    else
#	allWS="${allWS},${i}"; 
#    fi
#done 
#rundFtM  --ws ${allWS}