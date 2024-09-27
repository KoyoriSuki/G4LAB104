This program is used for gamma background evaluation under different height of lead/copper shielding

User defined UI commands:

/usrDefinedGeo/setCrystalType VALUE
	0: CMO 
	1: 4.5x4.5x4.5cm LMO 
	2: 1x1x1cm LMO
/usrDefinedGeo/setShieldingFlag VALUE
	0: without shielding
	1: shielding
/usrDefinedGeo/setShieldingHeight VALUE
	unit: cm
/usrDefinedGeo/setShieldingHeight2 VALUE
	unit: cm (higher shielding can have less thickness)
note: use /run/reinitializeGeometry after /usrDefinedGeo/ to correctly set the geometry

/mydet/setEnergyDistributionFlag VALUE
	not in use
/mydet/setSourceType VALUE
	0: use macro file settings
	1: randomly in shielding
	2: CRY generator

/histo/setFileName VALUE
	results will be output in ./VALUE_results/

Provided macro files in macFiles:
	init_vis.mac: for visualization. If run ./DBDecay without any parameters
	gammaBkg.mac: not in use
	shieldingX.mac: shielding with X cm in height
	na22.mac: counting rate simulation with Na-22 ion
	job.cmd: for condor submission
	mergeResults.cc: merge the result root files

To get the result:
	mkdir build -> cd build -> cmake .. -> make

	if run with job.cmd, the results will be:
		shieldingHeightResultX.root, X is the thread number
		copy the mergeResults.cc to ./build/resultFiles, then run "root mergeResults.cc", 
		it will give you the merged results

	or directly run "./DBDecay (xxx.mac) (X)"

To read the result:
	Examples in readResults.cc