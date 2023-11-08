function [R, X] = rQstar(Q, Qvop, compression, X)
% This code was written by Vincent Gras.
%
% Gras V, Boulant N, Luong M, Morel L, Le Touz N, Adam JP, Joly JC. 
% A mathematical analysis of clustering-free local SAR compression algorithms for MRI safety in parallel transmission. 
% IEEE transactions on medical imaging 2023;PP.
%
% Downloaded and slightly modified from:
% https://github.com/VincentGras/VOPcompressionCO
%
% Modification by Stephan Orzada: disable output to command window.
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

        X0 = sqrt( 1.0 / (SAR(Qvop, X0) + eps)) * X0;

        opt = optimset('Algorithm', 'sqp', 'MaxIter', 1000, ...
            'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on');

        if (compression)
            opt.OutputFcn = @stopCriterion;
        end

        X(:,i) = fmincon(@(x) objfun(Qi,x), X0, [], [], [], [], [], [], @(x) constrfun(Qvop, x), opt);
        
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
C = -V' * Q * V;
G = -2 * Q * V;

function [C,Ceq,G,Geq] = constrfun(Q,V)

Ceq = [];
C = reshape(pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none'), size(Q, 3), 1) - 1;
Geq = [];
G = 2 * reshape(pagemtimes(Q, V), numel(V), size(Q, 3));

function H = HessianFcn(Q, Qvop, V, lambda)

H = -2 * Q;

for c = 1:size(Qvop, 3)
    H = H + 2 * lambda.ineqnonlin(c) * Qvop(:,:,c);
end

function stop = stopCriterion(x, optimValues, state)

stop = (optimValues.fval < -1) && (optimValues.constrviolation <= 0);




