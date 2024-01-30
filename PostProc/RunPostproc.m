% This is an example on how to start the postprocessing.
% Written by Stephan Orzada at the German Cancer Research Center (DKFZ)

filename_in='/home/ubuntu/Data/SARMatrix_6x4_matrices_tri.mat'; %Points to the file containing the full set of SAR matrices.
name_VOP='VOP_Postproc_test_24ch_0p2'; %Points to a VOP-file (Contains VOPs and Sglobal).


load(name_VOP);
VOP_start=VOP;

[VOP,Sglobal_new]=VOP_overestimation_postprocessing_hybrid(filename_in,VOP_start,Sglobal,name_VOP);
