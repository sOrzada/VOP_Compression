function [r, x, y, out] = rco_fw_gt(Sv, A, compression , tol, maxit, y0)
% This code solves the CO-criterion (see Gras et al.) for SAR compression
% using Frank-Wolfe algorithm. It was developed using GLM 5.3 by Stephan
% Orzada @ German Cancer Research Center (DKFZ)
% This function can calculate r. It does not
% calculate the actual r if "compression" is true! In this case
% output values are 0 if below 1, 10 if above one
% and 5 if undecided.
% Sv: matrix under test
% A: set of VOPs as a cell array
% compression: if true, use break criterion. If false, calculate r to the
% end.
% tol: the tolerance for convergence. In practice 1e-3 to 1e-4 is sufficient
% maxit: the maximum number of iteration. I recommend > 500
% y0: can be used for a warm start. Be aware that numel(y0) is number of
% VOPs. This means that warm start is hardly feasible for compression, as
% the number of VOPs and therefore the number of elements in y grows.

if nargin<4 || isempty(tol),   tol=1e-3;  end %Tolerance of 1e-3 to 1e-4 is acceptable in practice and provides very fast results.
if nargin<5 || isempty(maxit), maxit=2000; end %Number of iterations is only crucial for very few matrices and can be kept low.
if nargin<6, y0=[]; end

options = optimset('TolX',1e-6,'Display','off');

thr=1; %Threshold for decision. >1 means that Sv is not dominated by VOPs
m = numel(A); %Number of VOPs
Nc = size(Sv,1); %Size of SAR matrices
Sv = (Sv+Sv')/2; %Ensure symmetry

A3 = zeros(Nc,Nc,m);

for w = 1:m
    A3(:,:,w) = (A{w}+A{w}')/2; 
end

A2 = reshape(A3, Nc*Nc, m);

%Prepare starting values for y
if isempty(y0)
    y = ones(m,1)/m;
else
    y = y0(:)/sum(y0); 
end

exceeds_thr=NaN; %This is the decision value, whether r is above one. NaN means undecided.
gap = inf; 
tau = 0; 
q = ones(m,1); 
x = zeros(Nc,1);

for k = 1:maxit
    Ay = reshape(A2*y, Nc, Nc);
    R  = chol(Ay);                       % Ay = R'*R
    B  = (R'\Sv)/R;
    B = (B+B')/2;       
    [U,E] = eig(B);  
    [tau,i] = max(real(diag(E)));
    x  = R\U(:,i);                       % x'Ay x = 1,  x'Sv x = tau
    q = real(pagemtimes(pagemtimes(x',A3),x));
    if compression
        if tau > thr*max(q),  exceeds_thr = true;  break; end   % r = tau/max(q) > thr: means definitely above 1
        if tau <= thr,        exceeds_thr = false; break; end   % r_CO <= tau <= thr: cannot be above 1
    end
    gap = tau*(max(q)-1);
    if gap <= tol*max(1,tau), break; end
    [~,ws] = max(q);
    a = fminbnd(@(al) fwF(y,ws,al,A2,Sv,Nc), 0, 1, ...
        options);
    y = (1-a)*y;  y(ws) = y(ws) + a;
end

x = x/sqrt(max(q));
if compression %If we use this function for compression, we are not interested in the actual value of r.
    if isnan(exceeds_thr) %Left loop without converging. Better safe than sorry, decide it is above 1.
        r=5; % Mark as undecided
    elseif exceeds_thr
        r=10; % Mark as over 1
    else
        r=0; %Mark as below 1
    end
else
    r = real(x'*Sv*x); % This is the actual value of r, which we provide if the function is not used for compression.
end
%disp(k)
out.ub = tau;  out.gap = gap;  out.y = y;
end

function t = fwF(y, ws, a, A2, Sv, Nc)
yy = (1-a)*y;  
yy(ws) = yy(ws)+a;
R  = chol(reshape(A2*yy, Nc, Nc));
B  = (R'\Sv)/R;  
B = (B+B')/2;         
t  = max(real(eig(B)));
end