function R = run_rco_fw(matrices, V_sub)
% This function is used to handle data between rco_fw_gt and the main
% compression algorithm function.

Nvop = size(V_sub,3);
Nmat=size(matrices,3);
Qvop_cell = cell(Nvop,1); %rco_fw_gt takes VOPs as cell array.
for a = 1:Nvop
    Qvop_cell{a}=(V_sub(:,:,a));
end

R=zeros(Nmat,1); %Initialize R

parfor a=1:Nmat %This parfor loop can be run in threaded environment
    [R(a), ~, ~, ~] = rco_fw_gt(matrices(:,:,a), Qvop_cell, true, 5e-4, 500);
end

% if any(R==5)
%     disp('Undecided')
% end


return
