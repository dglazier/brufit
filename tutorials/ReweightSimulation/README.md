# Reweight simulations so it mimics data more closely

Note your MC data must have an ID branch for this to work.
If it does not you may run,

    	root '$BRUFIT/macros/AddIDBranch.C("TreeName","FileName.root")'

which will have branch name UID (different from this example which generates
with fgID).


Run

	root CreateData.C
	
to produce a mock data and mock simulations file.

Run

	brufit CreateWeights.C
	
to create the weights file.

To use the simulation weights you append the following info
in the RooHSEventsPdf factory string :

"WEIGHTS@WeightSpeciesName,WeightFileName,WeightObjectName"

e.g. here it would be

"WEIGHTS@LikeData,Weights.root,MCWeights"