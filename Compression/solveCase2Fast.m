function [rCO, xopt, info] = solveCase2Fast(Sv, Bk, tol)
% Sv is a matrix to test
% Bk contains the VOPs
if nargin < 3, tol = 1e-9; end
K = numel(Bk);  Nc = size(Sv,1);

% ---- outer: projected gradient, BB steps, backtracking (ONE run) ----
beta = ones(K,1)/K;
[f,x,g] = evalOuter(beta,Sv,Bk);
fbest = f;  betabest = beta;  xbest = x;
bprev = beta;  gprev = g;  fprev = f;  eta = 1;      % <-- FIX: init fprev
for it = 1:400
    db = beta - bprev;  dg = g - gprev;  dd = real(db'*dg);
    if dd > 1e-14, eta = min(real(db'*db)/dd, 1e4); else, eta = 1; end
    bprev = beta;  gprev = g;
    for ls = 1:40                                   % backtracking
        btry = projSimplex(beta - eta*g);
        [ftry,xtry,gtry] = evalOuter(btry,Sv,Bk);
        if ftry <= f + 1e-4*real(g'*(btry-beta)) || eta < 1e-14, break; end
        eta = eta/2;
    end
    fprev = f;                                      % <-- FIX: store before update
    beta = btry;  f = ftry;  x = xtry;  g = gtry;
    if f < fbest, fbest = f;  betabest = beta;  xbest = x; end
    if it > 2 && abs(fprev - f) <= tol*max(1,abs(f)), break; end
end
UB = fbest;  beta = betabest;

% ---- top eigenspace at beta*: only if actually degenerate ----
Bbar = zeros(Nc); for k=1:K, Bbar = Bbar + beta(k)*Bk{k}; end
R = chol(Bbar);  C = (R'\Sv)/R;  C = (C+C')/2;
[W,L] = eig(C);  e = real(diag(L));  lam = e(end);
if lam - e(end-1) > 1e-8*max(1,abs(lam))
    xopt = xbest;                                   % simple eigenvalue: done
else                                                % targeted, vectorized
    Q  = orth(R \ W(:, e >= lam - 1e-8*max(1,abs(lam))));
    Cc = randn(size(Q,2),2000) + 1i*randn(size(Q,2),2000);
    Cc = [eye(size(Q,2)), Cc./sqrt(sum(abs(Cc).^2,1))];
    Xq = Q*Cc;  M = size(Xq,2);
    Sx = real(sum(conj(Xq).*(Sv*Xq),1));
    D  = zeros(K,M);
    for k = 1:K, D(k,:) = real(sum(conj(Xq).*(Bk{k}*Xq),1)); end
    [~,ix] = max(Sx ./ max(D,[],1));  xopt = Xq(:,ix);
end

% ---- achieved value + certificate ----
xopt = xopt/norm(xopt);
rCO  = real(xopt'*Sv*xopt) / max(cellfun(@(B) real(xopt'*B*xopt), Bk));
info.UB = UB;  info.beta = beta;
info.certified = (UB - rCO) <= 1e-6*max(1,abs(UB));
end

function [f,x,g] = evalOuter(beta,Sv,Bk)
K = numel(Bk);  Nc = size(Sv,1);
Bbar = zeros(Nc); for k=1:K, Bbar = Bbar + beta(k)*Bk{k}; end
R = chol(Bbar);  C = (R'\Sv)/R;  C = (C+C')/2;
[W,L] = eig(C);  e = real(diag(L));  [~,im] = max(e);
x = R\W(:,im);
s = real(x'*Sv*x);  d = real(x'*Bbar*x);
bk = zeros(K,1); for k=1:K, bk(k) = real(x'*Bk{k}*x); end
f = e(im);  g = -(s/d^2)*bk;                        % Danskin gradient
end

function p = projSimplex(v)
[u,~] = sort(v,'descend');  css = cumsum(u);  K = numel(v);
rho = find(u - (css-1)./(1:K)' > 0, 1, 'last');
p = max(v - (css(rho)-1)/rho, 0);
end