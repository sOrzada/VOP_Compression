% Example Code to use Hybrid Compression Algorithm
% Written by Stephan Orzada at the German Cancer Research Center (DKFZ)

%%%% User Defined Variables %%%%%
filename_in='/home/ubuntu/Data/SARMatrix_4x2.mat'; %Filename of matrices.

name_save = '4x2_20_sqrt0p5';


start_Overestimation = 20;  % Overestimation for first iteration in (%)
divider_Step = (0.5)^(1/2);         % After each iteration, the overstimation is multiplied with this factor
max_number_VOPs =2500;       % Maximum number of VOPs before program stops
max_iter = 6;               % Maximum number of iterations of the algorithm before it stops

%%%% Program %%%%

start_Overestimation = start_Overestimation /100;

VOP_compression_iterative_Hybrid(filename_in,start_Overestimation,divider_Step,[],name_save,2,max_number_VOPs,max_iter)