% Example Code to use Hybrid Compression Algorithm
% Written by Stephan Orzada at the German Cancer Research Center (DKFZ)

% This is an example of how to start the compression algorithm.
% There are two example datasets one in square and one in triangular
% format. (Square matrices can be reformatted to triangular form with the
% "format_sq2tri" function, either with single matrices or arrays of
% matrices, like matlab's "page*" functions.)
% The output of the compression function is a file containing the VOPs as
% well as other information on the compression, for example timings and the
% the positions of the VOPs in the original full set.

%%%% User Defined Variables %%%%%
filename_in='V_full.mat'; %Filename of matrices in square format.
%filename_in='V_full_tri.mat'; %Filename of matrices in triangular format.

start_Overestimation = 40;   % Overestimation for first iteration in (%) of SAR_wc
divider_Step = (0.5)^(1/3);  % After each iteration, the overstimation is multiplied with this factor

% Example for options:
options=[];
options.name_save='Example_8ch';
options.max_number_VOPs=160;
options.max_iter=4;
% options.continueFile='VOP_SOR_Example_8ch_0.2.mat';

% % The following are standard values, to which the function reverts, if no
% % other values are provided:

% options.useGPU=false;
% options.block_size_max=10000;
% options.block_size_2=50000;
% options.N_vop_switch=30;
% options.numMat4Pageeig=30000;
% options.numMat4GPU=100000;
% options.OverestimationMatrixType='Diagonal';
% options.OverestimationMatrix=[];




%%%% Program %%%%

start_Overestimation = start_Overestimation /100;

VOP_compression_iterative_Hybrid_tri(filename_in,start_Overestimation,divider_Step,options)