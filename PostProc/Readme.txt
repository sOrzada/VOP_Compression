New post processing algorithm for Virtual Observation Points (VOPs) written by Stephan Orzada at the German Cancer Research Center (DKFZ).

If you use this code, please cite: 
Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

This code combines the CC and the CO criteria for greatly enhanced speed.

When you read this, the paper describing this is not published, yet. In this case, please cite this as a combination of the methods from these two papers:

Orzada S, Fiedler TM, Quick HH, Ladd ME. Post-processing algorithms for specific absorption rate compression. Magn Reson Med 2021;86(5):2853-2861.

Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. IEEE Trans Med Imaging 2023;PP.

RunPostprocMulti.m is an example on how to call the Post Processing Algorithm.

VOP_overestimation_postprocessing_hybrid.m is the actual post processing algorithm.

rQstar.m was written by Vincent Gras for above paper. It has been modified to suppress the output to the command window.

The other scripts are helper functions.

If you encounter any problems or have any questions, do not hesitate to contact stephan.orzada@dkfz.de
