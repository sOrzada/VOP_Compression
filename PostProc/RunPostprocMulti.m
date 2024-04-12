% This script is used to start Post Processing of VOP files
% it is possible to load multiple files for postprocessing, but these must
% be from the same set of SAR matrices (same model).

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
    name_VOP=['PP_' name];
    VOP_start=VOP;

    [VOP,Sglobal_new]=VOP_overestimation_postprocessing_GPU_hybrid_tri(filename_in,VOP_start,Sglobal,name_VOP);
end