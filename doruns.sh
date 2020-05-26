#!/bin/bash


for run in 045 046 047 048 086 087 124 125 126 127 236 237 238 239 240
do
./FindBoltLocations /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar_TopInjector_PD3/BarrelSurveyFar_TopInjector_median/${run}.jpg /bigdisk/jamieson/TOW-Feb2020/BarrelSurveyFar_TopInjector_PD3/BarrelSurveyFar_TopInjector_median_texts/${run}.txt >& ${run}.log 
done



