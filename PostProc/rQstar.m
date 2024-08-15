function [R, X] = rQstar(Q, Qvop, compression, X)
% This code was originally written by Vincent Gras.
%
% Changes were introduced by Stephan Orzada:
% Reduced the number of unkowns by holding the phase of one channel
% constant. Reduced the tolerances for the optimizer and increased the
% number of function evaluations and iterations. Reduced the number of
% operations in Objective Function and Constraint Function. 
%
% Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. 
% A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. 
% IEEE transactions on medical imaging 2023;PP.
%
%
%
% maximizes xHQx subject to xHQ*x <= 1 
% Q           : SAR matrices
% Qvop        :  VOPs
% compression : Boolean value 
%               if true (e.g. when called by computeVOP_CO), then
%               optimization is stop as soon as xHQx exceeds 1
% X           : initial guess for the optimization
%               default is []
                  

if (nargin < 4)
    X = [];
end

if (nargin < 3)
    compression = false;
end

Nc = size(Qvop, 1);
N = size(Q, 3);

if (~isreal(Q) || ~isreal(Qvop))
    
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    [R, xr] = rQstar(CHtoRS(Q), CHtoRS(Qvop), compression, [real(X); imag(X)]);
    X = xr(1:Nc,:) + 1i * xr(Nc+1:2*Nc, :);
    
else
    
    if (isempty(X))
        X = zeros(Nc, N);
    end
    
    R = zeros(1, N);
    
    parfor i = 1:N
        
        Qi = double(Q(:,:,i));
        Qi = 0.5*(Qi+Qi');
        
        X0 = X(:, i);

        if (all(X0 == 0))
            [eigV, eigD] = eig(Qi, 'vector');
            eigD = real(eigD);
            [~, k] = max(eigD);
            X0 = real(eigV(:,k));

        end
        %IN the following lines we set the phase of the last channel to 0.
        Xt=X0(1:Nc/2,:)+1i*X0(Nc/2+1:Nc,:);
        Xt=Xt*exp(-1i*angle(Xt(end)));
        X0(1:Nc/2)=real(Xt);
        X0(Nc/2+1:Nc)=imag(Xt);

        X0 = sqrt( 1.0 / (SAR(Qvop, X0) + eps)) * X0;
        opt = optimoptions('fmincon',...
            'Algorithm', 'sqp', 'MaxIter', 2500, "MaxFunctionEvaluations",1500*Nc, ...
            'Display', 'off', 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, 'TolFun',1e-6, 'TolX',1e-6,...
            'StepTolerance',1e-12, 'ConstraintTolerance',1e-12, 'OptimalityTolerance', 1e-9);

        
        if (compression)
            opt.OutputFcn = @stopCriterion;
        end

        [X(:,i),~,~] = fmincon(@(x) objfun(Qi,x), X0, [], [], [], [], [], [], @(x) constrfun(Qvop, x), opt);


        R(i) = SAR(Qi,X(:,i))./SAR(Qvop,X(:,i));
        
        if (mod(i,1)== 0)
            %fprintf('parfor loop : %d / %d; %f\n', i, N, R(i));
        end
        
    end
    
end

function val = SAR(Q,V)
val = pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none');
val = max(val, [], 'all');
val = max(val, 0);

function [C, G] = objfun(Q, V)
V(end)=0; %Imaginary part of last channel can be kept 0.
temp=Q*V; %This can be used for Gradient AND Objective
C=-V'*temp;
G=-2*temp;
G(end)=0; %Gradient due to change in imaginary part of last channel has to be kept 0, so that the optimizer ignores this part.

function [C,Ceq,G,Geq] = constrfun(Q,V)
Ceq = [];
temp=pagemtimes(Q, V); %This can be used for Gradient and Constraint!
G=2*reshape(temp, numel(V), size(Q, 3));
C = reshape(pagemtimes(V, 'transpose', temp, 'none'), size(Q, 3), 1) - 1;
Geq = [];
G(end,:)=0;

function H = HessianFcn(Q, Qvop, V, lambda)

H = -2 * Q;

for c = 1:size(Qvop, 3)
    H = H + 2 * lambda.ineqnonlin(c) * Qvop(:,:,c);
end

function stop = stopCriterion(x, optimValues, state)

stop = (optimValues.fval < -1) && (optimValues.constrviolation <= 0);




