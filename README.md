# VOP_Compression

Speed improved algorithms for SAR compression.

If you use this code, please cite:
Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235

Needs Matlab 2023a or newer.

"Compression" contains a speed improved SAR compression algorithm that uses the CC[1] and CO[2] criteria in the iterative algorithm proposed by Orzada et al.[3]

"PostProc" contains a speed improved Virtual Observation Point (VOP) postprocessing algorithm that uses both the CC[1] and CO[2] critria within the post processing algorithm proposed by Orzada et al.[4]

Until the paper is published, please cite the appropriate original works by Orzada et al. and Gras et al.


1.	Lee J, Gebhardt M, Wald LL, Adalsteinsson E. "Local SAR in parallel transmission pulse design". Magn Reson Med 2012;67(6):1566-1578.
2.	Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. "A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission". IEEE Trans Med Imaging 2023;PP.
3.	Orzada S, Fiedler TM, Quick HH, Ladd ME. "Local SAR compression algorithm with improved compression, speed, and flexibility". Magn Reson Med 2021;86(1):561-568.
4.	Orzada S, Fiedler TM, Quick HH, Ladd ME. "Post-processing algorithms for specific absorption rate compression". Magn Reson Med 2021;86(5):2853-2861.


If you have any questions or suggestions, please contact stephan.orzada@dkfz.de
