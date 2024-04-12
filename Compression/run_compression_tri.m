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

name_save = 'Example_8ch';


start_Overestimation = 40;   % Overestimation for first iteration in (%)
divider_Step = (0.5)^(1/2);  % After each iteration, the overstimation is multiplied with this factor
max_number_VOPs =60;         % Maximum number of VOPs before program stops. The algorithm will change the step to target this exact number in the last iteration.
max_iter = 8;                % Maximum number of iterations of the algorithm before it stops

%%%% Program %%%%

start_Overestimation = start_Overestimation /100;

VOP_compression_iterative_Hybrid_double(filename_in,start_Overestimation,divider_Step,[],name_save,2,max_number_VOPs,max_iter)