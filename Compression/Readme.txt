Hybrid SAR matrix compression algorithm by Stephan Orzada from the German Cancer Research Center (DKFZ).

If you use this code, please cite: 
Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

The method is a combination of the CC criterion and CO criterion.

run_compression_tri.m is an example on how to use the code.

VOP_compression_iterative_Hybrid.m contains the actual compression algorithm. This can be run with or without a GPU. Check the function for more information.

rQstar.m was written by Vincent Gras for above paper. It has been modified to suppress the output to the command window.

The other scripts are helper functions.

If you encounter any problems or have any questions, do not hesitate to contact stephan.orzada@dkfz.de
