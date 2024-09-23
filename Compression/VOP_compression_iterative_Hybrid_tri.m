function VOP_compression_iterative_Hybrid_tri(filename_in,eps_G,R,options)
% 'filename_in' path to a .mat file with SAR matrices.
% 'eps_G'       Starting overestimation. 0.1 means 10% of SAR_wc
% 'R'           Reduction factor. After each iteration, overestimation is
%               multiplied by this factor.
% 'options'     Optional Struct with information for Algorithm:
%   'options.name_save'       Save name for results
%   'options.max_number_VOPs' Algorithm will try to exactly provide this
%                             number of VOPs.
%   'options.max_iter'        Maximum number of Iterations.
%   'options.OverestimationMatrixType'
%       'Diagonal' for a diagonal matrix. (Default)
%       'Global' for the mean of all SAR matrices. (~Global SAR).
%       'User' for a user defined matrix.
%   'options.OverestimationMatrix'  Provide and Overestimation Matrix,
%                                   either in NchxNch, or in triangular
%                                   format.
%   'options.continueFile'    filename of VOP file to continue.
%   'options.useGPU'          true: use a GPU, false: only use CPU.
%   'options.block_size_max'  Maximum number of matrices send to rQstar
%                             function in parallel. Default: 10,000.
%   'options.block_size_2'    Number of matrices send to rQstar function in
%                             parallel during 2nd step. Default: 50,000.
%   'options.N_vop_switch'    Number of VOPs at which criterion is switched
%                             from CC to CO. Default: 30.
%   'options.numMat4Pageeig'  Number of matrices send to pageeig function
%                             in parallel. Default: 50,000.
%   'options.numMat4GPU'      Number of matrices send to GPU. Use more if
%                             you have more Memory. Default: 100,000.
%
% This code was written by Stephan Orzada, German Cancer Research Center
% (DKFZ). 2023
% 
% This is a novel Hybrid approach with Lee's and Gras' criteria intermixed.
% Uses the iterative approach by Orzada et al., too.
%
% If you use this code, please cite: 
% Stephan Orzada, Thomas M. Fiedler, Mark E. Ladd. 
% "Hybrid algorithms for SAR matrix compression and the impact of post-processing on SAR calculation complexity", 
% Magn Reson Med 2024, https://doi.org/10.1002/mrm.30235
%
% filename_in: Points to a .mat file containting a variable "matrices" SAR matrices (Q-matrices); either (N_channel, N_channel, N_matrices) or lower triangular matrices (N_channel^2/2+N_channel/2, N_matrices).
% eps_G: Starting value for Overestimation, relative to worst case local SAR. (Note: The MINIMUM value of the overestimaton term is eps_G*SAR_wc)
% R: Factor to reduce the overestimation after each iteration step. (The overestimation factor eps_G is multiplied by R after each step)
% Sglobal: Provide Sglobal for algorithm (N_channel,N_channel). If you do not use it, provide anything and set "use_Sglobal" to 1 (for Sglobal according to Lee) or 2 for diagonal matrix.
% name_VOP: provide a string for a name for save-file.
% use_Sglobal: set to 0 to use your own Sglobal matrix, set to 1 to calculate it according to Lee, if you set this to 2 the algorithm will use a diagonal matrix.
% max_num_VOPs: maximum number of VOPs. This is the target for compression. The algorithm will reduce the overestimation adaptively to achieve exactly this number of VOPs, but will allow 5% more.
% 
% VOP contains the calculated VOPs (N_channel, N_channel, N_VOPs) including overestimation.
% The results are saved after each iteration in a file named "VOP_SOR_'name_VOP'_'eps_G'.mat". Here eps_G is the overestimation factor of respective iteration.
%
% SAR matrices and VOPs can be supplied in triangular form to save memory.
%
% This script needs R2023a or later.
% In this script you can choose whether you wish to do all calculations on
% CPU or have some calculations done on a GPU. (selected via the constant "use_GPU")
%
%  IMPORTANT NOTE:
%  The author takes no responsibility for the validity of the SAR models
%  generated with this tool. Validating the results is the responsibility
%  of the user!



%constants for internal purpose:

MaxFunEval=50000; %Maximum number of function evaluations in fmincon.
use_chol=false; %Use Cholesky Factorization if not run on GPU. If 'false' use pageeig if not run on GPU. Pageeig can be faster when many cores are available. I'm hoping there will be pagechol soon.

% In the following, we check the "options" input. If fields are not
% provided by the user, we revert to standard values.

if nargin<4
    warning('No options provided by user! Reverting to default values!')
end

try
    name_VOP=options.name_save;
catch
    name_VOP=filename_in;
end

try
    max_num_VOPs=options.max_number_VOPs;
catch
    warning('No maximum number of VOPs defined by user. Reverting to 100.')
    max_num_VOPs=100;
end

try
    max_iter=options.max_iter;
catch
    warning('No maximum number of Iterations defined by user. Reverting to 3.')
    max_iter=3;
end

try
    use_GPU=options.useGPU;
catch
    use_GPU=false; %Activate or deactivate the use of GPU. If set to "false" The algorithm will be run on CPU completely.
end

if use_GPU
    try
        gpuDevice
    catch
        warning('No useable GPU found. Reverting to pageeig on CPU.')
        use_GPU=false;
    end
end

try
    numMat4GPU=options.numMat4GPU;
catch
    numMat4GPU=100000; %Number of matrices that are transferred to GPU as a Batch
end

try
    numMat4Pageeig=options.numMat4Pageeig;
catch
    numMat4Pageeig=50000; %number of matrices used in the pageeig algorithm. Lower number reduces memory requirements.
end

try
    N_vop_switch=options.N_vop_switch;
catch
    N_vop_switch=100; %Number of VOPs at which we switch between Lee's criterion and Gras' criterion.
end

try
    block_size_max=options.block_size_max;
catch
    block_size_max=10000; %Number of matrices for rQstar calculation
end

try
    block_size_2=options.block_size_2;
catch
    block_size_2=50000; %Number of matrices for rQstar in second step.
end

try
    OverestimationMatrixType=options.OverestimationMatrixType;
catch
    OverestimationMatrixType='Diagonal';
end

if strcmp(OverestimationMatrixType,'User')
    try
        Sglobal=options.OverestimationMatrix;
        warning('No Sglobal provided. Reverting to Diagonal Matrix')
    catch
        OverestimationMatrixType='Diagonal';
    end
end

if strcmp(OverestimationMatrixType,'Global')
    use_Sglobal=1;
elseif strcmp(OverestimationMatrixType,'Diagonal')
    use_Sglobal=2;
elseif strcmp(OverestimationMatrixType,'User')
    use_Sglobal=3;
else
    error(['Unkown Matrix Type provided by user: ' OverestimationMatrixType ])
end


tic %Start timer for total elapsed time.

if R>=1
    error(['R must be < 1. Provided value is ' num2str(R)])
end
eps_decrease=R; %Had this as a constant initially. Was to lazy to rename all instances ;-)

if eps_G > 1
    eps_G=eps_G/100;
    warning(['epg_G is larger than 1. I assume this means you give it in %, so I divide by 100. New value is ' num2str(eps_G) '. Recommended starting values are 0.1 to 0.4.'])
end

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


if use_Sglobal==1 %If the user wants to use the actual Sglobal, calculate it
    Sglobal=sum(matrices,3);
    Sglobal=format_sq2tri(Sglobal);
elseif use_Sglobal==2 %If user wants a diagonal matrix, calculate it.
    Sglobal=format_sq2tri(diag(ones(y,1)));
elseif ismatrix(Sglobal)
    Sglobal=(Sglobal+Sglobal')/2; %make sure Sglobal is symmetric. (faster eig() and only real eigenvalues).
    Sglobal=format_sq2tri(Sglobal);
end

eigenvalues_max=zeros(N,1);
X=zeros(y,N); %Starting Vector for Vincent's optimization
disp('Calculating all Eigenvalues...')

parfor a=1:N %Using parfor is faster than using for in this case. If the user cannot use the parallelized form, "for" can be used instead. This is only a minor reduction in speed for large problems, since this part is only calculated ones.
    matrix_temp=matrices(:,a);
    for b=1:y
        matrix_temp(sq2tri(b,b))=real(matrix_temp(sq2tri(b,b))); %make sure that matrices are PSD.
    end
    matrices(:,a)=matrix_temp;
    current_mat=reformat_tri(matrix_temp,y);

    [V,eig_val_temp]=eig(current_mat,'vector');%Calculate maximum Eigenvalues
    
    [eigmax,Ieig]=max(real(eig_val_temp)); %Calculate maximum of real part.
    eigenvalues_max(a)=eigmax;
    X(:,a)=V(:,Ieig);
end

mean_matrix=mean(matrices,2);

disp('Sorting Eigenvalues...')
[B,I_eigen]=sort(eigenvalues_max,'descend'); %Sort maximum Eigenvalues. I_eigen contains the indexes.
I_eigen_save=I_eigen; %save I_eigen for later runs.
N_start=N; %Save number of matrices.


Sglobal=Sglobal/min(real(eig(reformat_tri(Sglobal,y))))*B(1); %Normalize Sglobal. We use the minimum eigenvalue here, because it determines the minimum value the overestimation term can assume, and this is crucial for the number of VOPs.
SAR_wc=B(1);
Sglobal=Sglobal*eps_G; %multiply with eps_G for overestimation control.
clear B %not needed anymore. We only need I, since it contains the order in which we look at the matrices Sv.


% This part checks, whether an existing VOP file is provided which will be
% continued.
try
    continueFile=options.continueFile;
    load(continueFile,'VOP','Sglobal','eps_G','VOPID')
    num_vops=size(VOP,3);
    V_sub=format_sq2tri(VOP);
    Sglobal=format_sq2tri(Sglobal);
    run_number=2; %Make this will lead to reduction of overestimation.
    iteration_step=1;
catch
    V_sub=matrices(:,I_eigen(1))+Sglobal; %Put first Voxel in subvolume and add Sglobal for Overestimation
    num_vops=1; %initialize num vops
    run_number=1;
    iteration_step=1;
    VOPID=I_eigen(1); %The index of the SAR-matrix with the highest eigenvalue is saved in VOPID
end

overestimation_factor=max(real(eig(reformat_tri(Sglobal,y))))/SAR_wc*100; %Calculate the maximum overestimation for information purposes.
overestimation_max=max(real(eig(reformat_tri(Sglobal,y))));
disp(['Starting: ' num2str(eps_G) ' Maximum Overestimation: ' num2str(overestimation_factor) '% of Worst Case SAR'])

too_many_vops=false; %to check whether a run has finished correctly or was broken due to too many vops.
ran_out_counter=0; %Counter for matrices which are included as VOPs without definitive descision.

options = optimoptions('fmincon','UseParallel',false,'Algorithm','sqp','display','off','SpecifyObjectiveGradient', true,'CheckGradients',false);
options.MaxFunctionEvaluations=MaxFunEval;

while (num_vops<max_num_VOPs) && (iteration_step<=max_iter)
    
    block_size=block_size_max;
    timings=[];
    inner_runs=0; %Number of runs in inner iteration
    I_eigen=I_eigen_save; %Since this is a new iteration step, we need all matrices again.
    I_eigen(1)=[]; %We know that the first matrix is in the VOPs, so we can delete it.
    num_vops_old=num_vops;
    V_sub=V_sub-Sglobal; %Due to memory and calculation time constraints, the sub set does not contain the original matrices as described in Lee's paper, but the matrices + overestimation. Here we recover the original matrices.
    eps_G_old=eps_G; %Save the value for target calculation
    if run_number>1 %If this is not the first iteration step, reduce the overestimation
        Sglobal=Sglobal*eps_decrease; %Change the overestimation term
        eps_G=eps_G*eps_decrease; %Change Overestimation factor (this is only used for display reasons hereafter).
    end
    disp(['Starting run #' num2str(run_number) '. eps_G = ' num2str(eps_G)])
    run_number=run_number+1;
    V_sub=V_sub+Sglobal; %Add the overestimation term to the matrices in the sub set. This saves calculation time, because it is only done once for each matrix in each iteration step.
    
    N=numel(I_eigen); %Total number of matrices in the full set.
    V_sub_square = reformat_tri(V_sub,y);
    while N>0 %Evaluate all matrices (voxels). 1st has already been put in.
        inner_runs=inner_runs+1;
        timings(1,inner_runs)=toc;
        timings(2,inner_runs)=N;
        a=1; %select next matrix.
        num_vops=size(V_sub,2); %Number of vops this round.

        if num_vops<=N_vop_switch || inner_runs==1
            
            Sv_current=matrices(:,I_eigen(a)); %Get current SAR matrix

            if num_vops>N_vop_switch
                disp('Start using CC criterion for first step...')
                Sv_current=mean_matrix;
            end

            A=zeros(num_vops,num_vops);
            A(1,:)=1; %this is to make sure that sum(c_wv)==1
            b=zeros(1,num_vops)';
            b(1)=1; %This is o make sure that sum(c_wv)==1
            c_wv_0=ones(1,num_vops)'/num_vops; %Starting values
            lb=zeros(num_vops,1); %Upper bound for c_wv
            ub=ones(num_vops,1); %lower bound for c_wv
            V_sub_double=double(V_sub); 
            Sv_current_double=double(Sv_current);

            [c_wv_out,feval,exitflag,~]=fmincon(@(c)optimize_me_tri(V_sub_double,Sv_current_double,c,y),c_wv_0,[],[],A,b,lb,ub,[],options);

            c_wv=reshape(c_wv_out,[1 num_vops]); %reshape c_wv for easier multiplication with matrices.

            if num_vops>N_vop_switch %If this is just a pre-test, we check mean_matrix, which is not in the actual SAR matrix set.
                feval=-1;
            end

            if feval>0 %If not upperbounded, add it to subset.
                if exitflag==0 %inform the user that a VOP was added that ran out of iterations. We only do this here, since it might cause unnecessary VOPs to be added. Stoping prematurely when the result is PSD will just lead to less dropped matrices in Chol-step, but has no negative effect on results.
                    ran_out_counter=ran_out_counter+1;
                    disp('Ran out of iterations. Consider increasing "MaxFunEval".')
                end
                if num_vops<=N_vop_switch
                    V_sub(:,num_vops+1)=Sv_current+Sglobal; %Add current matrix to subset.
                    disp([num2str(eps_G) ', Finished: ' num2str((N_start-N+a)/N_start*100) '% VOP: ' num2str(num_vops+1) '. ID: ' num2str(I_eigen(1))])
                    VOPID = [VOPID I_eigen(1)]; %Since the matrix with the index I_eigen(1) was found to be a VOP, place its ID in VOPID
                    I_eigen(1)=[];
                    num_vops=size(V_sub,2);
                end
                V_sub_square = reformat_tri(V_sub,y);
            else
                P_pre2=sum(double(V_sub).*(c_wv/sum(c_wv)),2); %Calculate the weighted sum of the VOPs here
                if num_vops<=N_vop_switch
                    I_eigen(1)=[];
                end
                %Pretest: Check main Diagonal of P_pre2-Sv for Elements<0.
                %These are not PSD and do not need to be checked with cholesky
                %decomposition. This check is extremely fast and saves a lot of
                %time especially for large channel counts.
                x=1:y;
                pos=(x-1).*x./2+x; %position in vector of lower triangle
                N=numel(I_eigen);
                isPSD=true(N,1);
                isPSD=isPSD & not(any(P_pre2(pos)-matrices(pos,I_eigen)<0,1)).'; %All matrices with main diagonal <=0 are not PSD.
                disp([num2str((N-sum(double(isPSD)))/N*100) '% of matrices not PSD in pre-test']) %Just to check whether it does anything!
                I_eigen_unchecked_pos=find(isPSD); %Position in I_eigen. We need this later to delete unnecessary positions in I_eigen
                I_eigen_unchecked=I_eigen(isPSD);

                if use_GPU==true %GPU enabled version
                    disp('Check for PSDness using Cholesky Factorization on GPU...')
                    N_unchecked=numel(I_eigen_unchecked); %Now we only need to check those matrices we haven't excluded in the step before.
                    N_save=N;
                    gpu_Steps=[1:numMat4GPU:N_unchecked N_unchecked]; %Steps for Transmission to GPU
                    isPSD_complete=false(N_unchecked,1); %Allocate memory
                    for a=1:numel(gpu_Steps)-1
                        [isPSD]=check_PSD_chol_tri(gpuArray(P_pre2)-gpuArray(matrices(:,I_eigen_unchecked(gpu_Steps(a):gpu_Steps(a+1)))),y); %Check for PSDness on GPU
                        isPSD_complete(gpu_Steps(a):gpu_Steps(a+1))=isPSD; %Save the result chunk in total result
                    end
                    delete_pos=double((isPSD_complete)); %Just for userinformation.
                    dropped_mats=sum(delete_pos); %Just for user information.
                    I_eigen(I_eigen_unchecked_pos(isPSD_complete))=[]; %Delete the entries from this vector. Only matrices which this vector still points to are included in further calculations.
                    N=numel(I_eigen); %Determine new number of matrices we are still working with.
                    disp(['Dropped ' num2str(dropped_mats) ' of ' num2str(N_save) ' matrices (Cholesky).'])
                
                elseif use_chol %Version for CPU using Cholesky Factorization
                    disp('Check for PSDness using Cholesky Factorization on CPU...')
                    kappa=0;
                    N_save=N;
                    delete_pos=zeros(N,1);
                    while kappa<N %This is not a parfor loop since "chol" doesn't work well in a parallel loop.
                        kappa=kappa+1;
                        if isPSD(kappa)
                            Sv_current=matrices(:,I_eigen(kappa));
                            P=P_pre2-Sv_current;
                            [~,flag_chol] = chol(reformat_tri(P,y)); %The triangular matrix needs to be reformatted so that chol() can work with it.
                            if flag_chol==0
                                delete_pos(kappa)=1;
                            end
                        end
                    end
                    dropped_mats=sum(delete_pos); %This is for info only.
                    I_eigen(logical(delete_pos))=[];
                    N=numel(I_eigen);
                    disp(['Dropped ' num2str(dropped_mats) ' of ' num2str(N_save) ' matrices using coefficients from last step.'])
                else %Version for CPU using pageeig. Better CPU utilization than chol, especially with many cores. Hopefully Matlab will contain pagechol in a future version.
                    disp('Check for PSDness using pageeig...')
                    N_unchecked=numel(I_eigen_unchecked); %Now we only need to check those matrices we haven't excluded in the step before.
                    N_save=N;
                    cpu_Steps=[1:numMat4Pageeig:N_unchecked N_unchecked]; %Steps for blockwise calculation
                    isPSD_complete=false(N_unchecked,1); %Allocate memory
                    
                    for a=1:numel(cpu_Steps)-1
                        eigenv=pageeig(double(reformat_tri(P_pre2-double(matrices(:,I_eigen_unchecked(cpu_Steps(a):cpu_Steps(a+1)))),y)));
                        min_eig=squeeze(min(eigenv));
                        isPSD=min_eig>=0;
                        isPSD_complete(cpu_Steps(a):cpu_Steps(a+1))=isPSD; %Save the result chunk in total result
                    end
                    delete_pos=double((isPSD_complete)); %Just for user information.
                    dropped_mats=sum(delete_pos); %Just for user information.
                    I_eigen(I_eigen_unchecked_pos(isPSD_complete))=[]; %Delete the entries from this vector. Only matrices which this vector still points to are included in further calculations.
                    N=numel(I_eigen); %Determine new number of matrices we are still working with.
                    disp(['Dropped ' num2str(dropped_mats) ' of ' num2str(N_save) ' matrices using coefficients from last step.'])
                end
            end
        elseif inner_runs>2
            
            if numel(I_eigen)<block_size
                block_size_current = numel(I_eigen);
            else
                block_size_current = block_size;
            end
            [R, X_temp] = rQstar(reformat_tri(matrices(:,I_eigen(1:block_size_current))/SAR_wc,y), V_sub_square/SAR_wc, true, X(:,I_eigen(1:block_size_current)));
            X(:,I_eigen(1:block_size_current))=X_temp;

            if all(R<=1)
                I_eigen(1:block_size_current)=[];
                block_size=block_size_max;
            else
                I_eigen(R<=1)=[];
                drop_percent=sum(double(R<=1))/block_size;
                block_size=floor(block_size_max*drop_percent)+1; %Adaptive blocksize reduces multiple checks of same matrices.
                V_sub(:,num_vops+1)=matrices(:,I_eigen(1))+Sglobal;
                disp([num2str(eps_G) ', Finished: ' num2str((N_start-N+a)/N_start*100) '% VOP: ' num2str(num_vops+1) '. ID: ' num2str(I_eigen(1)) ' Drop: ' num2str(drop_percent*100) '%'])
                VOPID = [VOPID I_eigen(1)]; %Since the matrix with the index I_eigen(1) was found to be a VOP, place its ID in VOPID
                I_eigen(1)=[];
                V_sub_square = reformat_tri(V_sub,y);
                num_vops=size(V_sub,2);
            end
        else %In the second step, throw out all remaining dominated matrices, while we still have a low number of VOPs.
            disp('Start using CO criterion on all remainging matrices...')
            N=numel(I_eigen);
            check2_Steps=[1:block_size_2:N N];
            delete_pos=false(N,1);
            for a=1:numel(check2_Steps)-1
                [R, X_temp] = rQstar(reformat_tri(matrices(:,I_eigen(check2_Steps(a):check2_Steps(a+1)))/SAR_wc,y), V_sub_square/SAR_wc, true, X(:,I_eigen(check2_Steps(a):check2_Steps(a+1))));
                X(:,I_eigen(check2_Steps(a):check2_Steps(a+1)))=X_temp;
                delete_pos(check2_Steps(a):check2_Steps(a+1)) = R<=1;
            end
            dropped_mats=sum(double(delete_pos));
            I_eigen(delete_pos)=[];
            disp(['Dropped ' num2str(dropped_mats) ' matrices in 2nd step.'])
        end
        N=numel(I_eigen);
        if mod(inner_runs,1000)==0
            disp(['Finished: ' num2str((1-N/N_start)*100) '%'])
        end 
    end
    elapsed_time=toc;
    %Save the time since start and number of remaining VOPs for user
    %information
    timings(1,inner_runs)=toc;
    timings(2,inner_runs)=N;
    if not(too_many_vops)
        VOP=zeros(y,y,num_vops); %Allocate memory so it doens't have to grow in the loop.
        for x1=1:num_vops
            VOP(:,:,x1)=single(reformat_tri(V_sub(:,x1),y)); %place all matrices in VOP variable
        end
        Sglobal=single(reformat_tri(Sglobal,y));
        %disp([num2str(ran_out_counter/num_vops*100) '% VOPs with insufficient iterations for descision.'])
        save(['VOP_SOR_' name_VOP '_' num2str(eps_G) '.mat'],'VOP','Sglobal','elapsed_time','eps_G','overestimation_max','VOPID','timings') %save VOPs in file with a name that tells what was compressed and how.
        Sglobal=format_sq2tri(Sglobal);
    else
        disp('Too many VOPs. Finishing.')
        %save(['VOP_SOR_' name_VOP '_' num2str(eps_G) '.mat'],'VOP','Sglobal','elapsed_time','eps_G','ran_out_counter','overestimation_max') %save VOPs in file with a name that tells what was compressed and how.
    end
    % The following part adapts the reduction factor R to get as close as
    % possible to max_num_vops.
    num_vops=size(VOP,3);
    vop_factor=num_vops/num_vops_old; %Calculate the factor by which the VOPs have increased during latest iteration. This will be very nearly the same factor for the next iteration.
    estimated_number_of_VOPs=vop_factor*num_vops; %The number of VOPs the next iteration will bring, if R remains the same.
    if iteration_step>1 && estimated_number_of_VOPs > max_num_VOPs %The first iteration does not give a valid result.
        slope=(log10(eps_G)-log10(eps_G_old))/(log10(num_vops)-log10(num_vops_old)); %Empirically I found that in loglog-representation the overestimation over the number of VOPs is a straight line (for more than 10-20 VOPs), so we can approximate the behaviour.
        offset=log10(eps_G_old)-slope*(log10(num_vops_old));
        new_eps_g=10^(log10(max_num_VOPs)*slope+offset); %Calculate the new overestimation which will most likely result in the desired number of VOPs.
        eps_decrease=new_eps_g/eps_G;
        disp(['*** Changed R to ' num2str(eps_decrease) ' ***'])
    end
    iteration_step=iteration_step+1;
    
    disp(['Finished with ' num2str(size(VOP,3)) ' VOPs. Elapsed time is: ' num2str(uint32(elapsed_time)) ' s.'])
end

end %End function

function [C,G]=optimize_me_tri(B,Bi,c_wvb,ch) %This is the function for optimization by fmincon
%This is the function for optimization by fmincon. This uses only lower
%triangles to increase speed and reduce memory requirements
[V,eig_values]=(eig(reformat_tri(Bi-sum(B.*c_wvb',2),ch),'vector')); %calculate eigenvalues
[C,I]=max(eig_values); %We want to minimize the maximum eigenvalue. Ideally, all Eigenvalues are negative, so that Bi is dominated.
G=-real(pagemtimes(V(:,I),'ctranspose',pagemtimes(reformat_tri(B,ch),V(:,I)),'none'));

end

