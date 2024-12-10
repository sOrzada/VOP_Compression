New post processing algorithm for Virtual Observation Points (VOPs) written by Stephan Orzada at the German Cancer Research Center (DKFZ).

If you use this code, please cite: 
Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

This code combines the CC and the CO criteria for greatly enhanced speed.

RunPostprocMulti.m is an example on how to call the Post Processing Algorithm.

VOP_overestimation_postprocessing_hybrid.m is the actual post processing algorithm.

rQstar.m was written by Vincent Gras for above paper. It has been modified to suppress the output to the command window.

The other scripts are helper functions.

If you encounter any problems or have any questions, do not hesitate to contact stephan.orzada@dkfz.de
