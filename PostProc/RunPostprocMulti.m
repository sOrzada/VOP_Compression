% This script is used to start Post Processing of VOP files
% it is possible to load multiple files for postprocessing, but these must
% be from the same set of SAR matrices (same model).

% Example for options for post processing algorithm:
options=[];
options.use_GPU=false;
options.N_block=400;
options.numMat4GPU=100000;
options.numMat4Pageeig=10000;
options.sort_order='ascend';
% options.name is set in the loop below.
% There are more options. Please refer to the function for more
% information.

[file, path]=uigetfile('.mat','Select SAR matrices','MultiSelect','off');

filename_in=[path '/' file];

[file, path]=uigetfile('.mat','Select Compression Results','MultiSelect','on');

if iscell(file)
    num_files=numel(file);
else
    num_files=1;
end

for a=1:num_files
    
    if iscell(file)
        load([path '/' file{a}])
        filename=file{a};
    else
        load([path '/' file])
        filename=file;
    end
    [~,name,~]=fileparts(filename);
    Nch=size(Sglobal);
    options.name_save=name;
    VOP_start=VOP;

    [VOP,Sglobal_new]=VOP_overestimation_postprocessing_GPU_hybrid_tri(filename_in,VOP_start,Sglobal,options);
end