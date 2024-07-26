Hybrid SAR matrix compression algorithm by Stephan Orzada from the German Cancer Research Center (DKFZ).

If you use this code, please cite: 
Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

The method is a combination of the CC criterion and CO criterion.
When you read this, the paper describing this is not published, yet. In this case, please cite this as a combination of the methods in these two papers:

Orzada S, Fiedler TM, Quick HH, Ladd ME. Local SAR compression algorithm with improved compression, speed, and flexibility. Magn Reson Med 2021;86(1):561-568.

Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. IEEE Trans Med Imaging 2023;PP.

run_compression_tri.m is an example on how to use the code.

VOP_compression_iterative_Hybrid.m contains the actual compression algorithm. This can be run with or without a GPU. Check the function for more information.

rQstar.m was written by Vincent Gras for above paper. It has been modified to suppress the output to the command window.

The other scripts are helper functions.

If you encounter any problems or have any questions, do not hesitate to contact stephan.orzada@dkfz.de
