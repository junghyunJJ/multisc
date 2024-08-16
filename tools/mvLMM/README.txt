
mvLMM - matrix-variate linear mixed-model is a linear mixed model-solver for multiple phenotypes.
It requires the package pylmm (https://github.com/nickFurlotte/pylmm).
pylmm can be used to read in PLINK formatted files.
Kinship matrices can be calculated with any available program such as GCTA or fastLMM but will need to be converted to the correct format.  Alternatively, the user can calculate their own kinship matrix.

With a matrix of phenotypes (nx2 for example) and an nxn kinship matrix, you simply create an mvLMM object, perform the optimization and then you have the results.

Check out example.py.
Please email nick.furlotte@gmail.com with any questions.
