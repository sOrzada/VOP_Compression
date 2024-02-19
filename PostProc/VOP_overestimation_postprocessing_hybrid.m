function [VOP,Sglobal_new]=VOP_overestimation_postprocessing_hybrid(filename_in,VOP_start,Sglobal,name_VOP)
% This code was written by Stephan Orzada, German Cancer Research Center (DKFZ).
% In 2023; email: stephan.orzada@dkfz.de
% This code reduces overestimation for a given set of VOPs.
%
% This is a new version, enhancing the original post processing by applying
% the CO criterion by Vincent Gras et al.
% Paper is in preperation. Until publication please cite this work a 
% combination of the following two papers:
%
% Orzada S, Fiedler TM, Quick HH, Ladd ME. 
% Post-processing algorithms for specific absorption rate compression. 
% Magn Reson Med 2021;86(5):2853-2861.
%
% Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. 
% A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. 
% IEEE Trans Med Imaging 2023;PP.
%
% 'filename_in' points to a file containing the full set of SAR-matrices
% (named 'matrices' and is [Nchannel,Nchannel,Nmatrices])
% 'VOP_start' are the VOPs from an earlier compression (Using Lee's or
% Orzada's compression algorithms). It is assumed here that these contain
% overestimations.
% 'Sglobal' is the orignal overestimation matrix used in compression of the
% VOPs
% 'name_VOP' is a string containing a name for the save-file
%
% When using a GPU, it is recommend to use matrices in single precision.
% 
% This code uses a hyrbid approach with the "rQstar" calculation by V.Gras et al.
%
%  IMPORTANT NOTE:
%  The author takes no responsibility for the validity of the SAR models
%  generated with this tool. Validating the results is the responsibility
%  of the user!


%%%% TODO: Use Chol in parfor loop (is this really neccessary?)

tic %Start timer for total Elapsed time.
%constants for internal purpose:
R=1; % (0 to 1) how much of the original overestimation do we subtract from the VOPs? (Experimental. Better leave at 1)
exponent_cwv=1; %can change weighting of coefficients for adding overestimation. (Experimental Better leave at 1)
use_global_opt=false; %Use global optimization. Very slow. Only minimal enhancement if any.
use_parallel=false; %Not sure yet whether this makes it any faster. (Might be system dependent)
use_outputfun=false; %Stops optimization when PSD. Not a good idea, since fewer matrices are dropped in cholesky step.
use_GPU=false; %If you have a CUDA GPU which can fit your matrices in its memory, try doing the cholesky factorization on the GPU. I did it on a GV100 with 32GB and it was incredibly fast.
N_block=400; %Number of Matrices that are checked in parallel with Vincent's criterion.
numMat4GPU=2000000; %Number of matrices that are transferred to GPU as a Batch
sort_order='ascend'; %Order in which matrices are sorted. Use 'ascend' or 'descend'. 'ascend' seems to provide better results.

disp('Loading matrices...')
load(filename_in,'matrices');
disp('Loading complete.')
try
    dim_count=ndims(matrices);
    if dim_count==3
        matrices=format_sq2tri(matrices);
    end
    [x,N]=size(matrices); %Determine size of Omega
    y=-0.5+sqrt(0.5^2+2*x);%Calculate number of channels from number of elements in first dimension
catch
    error('File does not contain a variable named "matrices"!')
end

try
    if size(Sglobal,2)>1
        Sglobal=format_sq2tri(Sglobal);
    end
catch
    error('Please Provide Sglobal in a correct way!')
end

if ndims(VOP_start)==3
    VOP_start=format_sq2tri(VOP_start);
end

eigenvalues_max=zeros(N,1);

disp('Calculating all Eigenvalues...')

parfor a=1:N %Using parfor is faster than using for in this case. If the user cannot use the parallelized form, "for" can be used instead. This is only a minor reduction in speed for large problems, since this part is only calculated ones.
    matrix_temp=matrices(:,a);
    for b=1:y
        matrix_temp(sq2tri(b,b))=real(matrix_temp(sq2tri(b,b))); %make sure that matrices are PSD.
    end
    matrices(:,a)=matrix_temp;
    current_mat=reformat_tri(matrix_temp,y);
    eig_val_temp=eig(current_mat);%Calculate maximum Eigenvalues
    eigenvalues_max(a)=max(real(eig_val_temp)); %Calculate maximum of real part.
end

disp('Sorting Eigenvalues...')
[B,I_eigen]=sort(eigenvalues_max,sort_order); %Sort maximum Eigenvalues. I_eigen contains the indexes.
%[~,I_eigen_resort]=sort(I_eigen,'ascend'); 

disp('Reordering matrices...')

if strcmp(sort_order,'ascend')
    max_eig=B(I_eigen(end));
else
    max_eig=B(I_eigen(1));
end

num_vops=size(VOP_start,2); %initialize num vops


clear B %not needed anymore. We only need I, since it contains the order in which we look at the matrices Sv.


V_sub=VOP_start-Sglobal*R; %Take VOPs and reduce to original matrices.
V_sub_pure=V_sub;   %Save the VOPs without overestimation in this variable, we will need it later.
Sglobal_new=V_sub*0+Sglobal*(1-R); %Allocate individual overestimation matrices.

mean_matrix=mean(matrices,2); %Mean of all matrices. This will elt us drop more matrices in first step.

N_start=N; %Save number of matrices.
a=0;
runs=0;


%Prepare constraints
A=zeros(num_vops,num_vops);
A(1,:)=1; %Together with b(1)=1 this makes sure that sum(c_wv)==1
b=zeros(1,num_vops)';
b(1)=1;
c_wv_0=ones(1,num_vops)'/num_vops; %Startvalue for optimization
lb=zeros(num_vops,1); %Lower bound for c_wv
ub=ones(num_vops,1);  %upper bound for c_wv

while a<N %Evaluate all matrices (voxels).
    a=1;
    runs=runs+1;
    Sv_current=matrices(:,I_eigen(a)); %Get current SAR matrix
 

    V_sub_double=double(V_sub); %fmincon needs double values. Here we use the sub volume WITH overestimation to see whether this is sufficient to dominate the matrix we are looking at.

    
    if runs==1 %First run using Lee's criterion with Kuehne's speed-up. Throws out ~50% of the matrices very quickly.
        
       
        disp('Calculating coefficients for first drop.')
        if use_outputfun==true
            options = optimoptions('fmincon','UseParallel',use_parallel,'OutputFcn',@outfun,'Algorithm','sqp','display','off'); %Turn of display. Otherwise notifications will be overwhelming
        else
            options = optimoptions('fmincon','UseParallel',use_parallel,'Algorithm','sqp','display','off'); %Turn of display. Otherwise notifications will be overwhelming
        end
        Sv_current_double=double(mean_matrix);

        [c_wv_out,~]=fmincon(@(c)optimize_me_tri(V_sub_double,Sv_current_double,c,y),c_wv_0,[],[],A,b,lb,ub,[],options);
        c_wv=reshape(c_wv_out,[1 num_vops]); %reshape c_wv for easier multiplication with matrices.
        
        %Here, we do not care whether this matrix is upperbounded. We just
        %want the coefficients to through out as many matrices as possible.
        can_be_upperbounded=1; %...set to 1.
        

    else %All further runs are performed on just a block of matrices, and Overestimation is calculated when matrices are not already dominated. Many matrices with a low index in I_eigen can be dominated!
        
        if N<N_block
            N_block=N;
        end
        [R, ~] = rQstar(reformat_tri(matrices(:,I_eigen(1:N_block)),y), reformat_tri(V_sub_double,y), true,[]);
        I_eigen(R<1)=[];
        N=numel(I_eigen);
        
        if all(R<1)
            can_be_upperbounded=1;
        else
            can_be_upperbounded=0;

        end
    end
    if N>0
        Sv_current=matrices(:,I_eigen(1));
    end
    
    %Calculate coefficients for adding overestimation.
    if can_be_upperbounded==0 %cannot be dominated. We need to add some more overestimation.
        %In the first step, we find the coefficients c_wv that minimize the
        %"distance" between Sv_current and the "pure" VOPs (without
        %overestimation).

        V_sub_double=double(V_sub_pure); %Here we use the sub volume WITHOUT overestimation, so that we get a result with minimum total overstimation
        Sv_current_double=double(Sv_current);
        if use_global_opt==0 %fmincon with one start value
            options = optimoptions('fmincon','UseParallel',use_parallel,'Algorithm','sqp','display','off');
            [c_wv_out,~]=fmincon(@(c)optimize_me_tri(V_sub_double,Sv_current_double,c,y),c_wv_0,[],[],A,b,lb,ub,[],options);
            c_wv=reshape(c_wv_out,[1 num_vops]); %reshape vor easier multiplication with matrices.
        end
        
        %In the next step, we calculate what we have to add to the VOPs in
        %order to dominate Sv_current.
        c_wv=c_wv/sum(c_wv); %normalue c_wv (might be unnecessary, but doesn't cost us much.)
        %Update Sglobal, so that Sv_current is dominated.
        for axx=1:1 %Could be run several times to avoid numerical error, but no problems have occured during my tests.
            P1_pre=sum(V_sub_pure.*(c_wv/sum(c_wv)),2); %weight sub volume with coefficients.
            [V,E]=eig(reformat_tri(P1_pre-Sv_current,y));
            E=real(E); %E should be real. It mostly isn't due to numerical error
            E_plus=E;
            E_plus((E<0))=0; %Caculate E+
            E_minus=E_plus-E; %Calculate E-
            Z=V*E_minus*V'; %Calculate the Z-matrix so that the matrix we are looking at is dominated
            Z=format_sq2tri((Z+Z')/2); %Make sure Z is symmetric (due to numerical error it mostly isn't.)
            for kx=1:num_vops %calculate the additonal overestimation for each VOP.
                [V,E]=eig(reformat_tri(Sglobal_new(:,kx)-Z*((c_wv(:,kx)^exponent_cwv)/sum(c_wv.^(exponent_cwv+1))),y)); %here we weight each VOP according to its coefficient.
                E=real(E); %E should be real. It mostly isn't due to numerical error
                E_plus=E;
                E_plus((E<0))=0; %Caculate E+
                E_minus=E_plus-E; %Calculate E-
                Z2=V*E_minus*V'; %Calculate the Z-matrix so that the matrix we are looking at is dominated
                Z2=(Z2+Z2')/2; %Make sure Z is symmetric (due to numerical error it mostly isn't.)
                Sglobal_new(:,kx)=Sglobal_new(:,kx)+format_sq2tri(Z2); %Update S_global_new
                V_sub(:,kx)=V_sub_pure(:,kx)+Sglobal_new(:,kx); %Update sub volume.

            end

        end
        %Calculate maximum value of overestimation term. (The user will
        %want to know)
        max_Value_S=0;
        I_eigen(1)=[]; %We have integrated enough overestimation to include this matrix. We do not need to check it again.
        N=numel(I_eigen);
        for VOP_counter=1:num_vops
            max_Value_S_temp=max(eig(reformat_tri(Sglobal_new(:,VOP_counter),y)));
            if max_Value_S_temp>max_Value_S
                max_Value_S=max_Value_S_temp;
            end
        end
        eps_G=max_Value_S/max_eig;
        disp([num2str(eps_G) ', Finished: ' num2str((N_start-N+a)/N_start*100) '%.'])
    end
    
    %Check for matrices that are dominated by the previous result. Saves a
    %lot of calculateion time, because we do not have to do the fmincon
    %optimization again.
    P_pre2=sum(V_sub.*(c_wv/sum(c_wv)),2);
    
    if use_GPU==true && runs==1  %Do the Cholesky factorization on a GPU
        disp('Start dropping matrices with Cholesky Factorization on GPU.')
        

        N_save=N;

        gpu_Steps=[1:numMat4GPU:N N];
        isPSD_complete=false(N,1);
        for kappa=1:numel(gpu_Steps)-1
            gMatrices=gpuArray(P_pre2)-gpuArray(matrices(:,I_eigen(gpu_Steps(kappa):gpu_Steps(kappa+1))));
            [isPSD]=check_PSD_chol_tri(gMatrices,y); %Here we use a massively parallel implementation of cholesky on a GPU.
            isPSD_complete(gpu_Steps(kappa):gpu_Steps(kappa+1))=isPSD;
        end
        delete_pos=double((isPSD_complete));
        dropped_mats=sum(delete_pos); %Just for user information.
        I_eigen(isPSD_complete)=[]; %Delete the entries from this vector. Only matrices which this vector still points to are included in further calculations.
        N=numel(I_eigen); %Determine new number of matrices we are still working with.
        disp(['Dropped ' num2str(dropped_mats) ' of ' num2str(N_save-a) ' matrices using coefficients from last step.'])
    elseif runs==1 %Do the Cholesky factorization on the CPU.
        disp('Start dropping matrices with Cholesky Factorization on CPU.')
        
        N_save=N;
        delete_pos=zeros(N,1);
        
        for kappa=1:N %Parfor wouldn't speed up much.

            Sv_current=matrices(:,I_eigen(kappa));
            P=P_pre2-Sv_current;
            [~,flag] = chol(reformat_tri(P,y));
            if flag==0
                delete_pos(kappa)=1;
            end
        end


        dropped_mats=sum(delete_pos);
        I_eigen(logical(delete_pos))=[];
        N=numel(I_eigen);
        disp(['Dropped ' num2str(dropped_mats) ' of ' num2str(N_save-a) ' matrices using coefficients from last step.'])
        
    end
end
elapsed_time=toc;
disp(['Finished with ' num2str(num_vops) ' VOPs. Elapsed time is: ' num2str(uint32(elapsed_time)) ' s.'])

VOP=reformat_tri(V_sub,y); %place all matrices in VOP variable
Sglobal=reformat_tri(Sglobal,y);
Sglobal_new=reformat_tri(Sglobal_new,y);

save(['VOP_PP_SOR_' name_VOP '_' num2str(eps_G) '.mat'],'VOP','Sglobal','Sglobal_new','elapsed_time','eps_G','max_Value_S') %save VOPs in file with a name that tells what was compressed and how.
end

function out=optimize_me_tri(B,Bi,c_wvb,ch) %This is the function for optimization by fmincon
%This is the function for optimization by fmincon. This uses only lower
%triangles to increase speed and reduce memory requirements
eig_values=(eig(reformat_tri(sum(B.*c_wvb',2)-Bi,ch))); %calculate eigenvalues
out=-min(eig_values); %negativ of minimum eigenvalue (fmincon minimizes, but we want to maximize the minimum eigenvalue).
end

function stop = outfun(~,optimValues,~) %Stop when you see it is dominated. Should not do this, since the number of VOPs dropped due to cholesky shrinks.
    if optimValues.fval < 0
        stop = true;
    else
        stop = false;
    end
end

