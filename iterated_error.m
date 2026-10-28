function [err,vals] = iterated_error(sol,x,fRef,kappaFun,alphaFun,betaFun,nq)
%ITERATED_ERROR  Sup-norm error of the iterated approximation.
%
%   The function returns the error of the iterated collocation or iterated
%   modified solution
%
%       ftilde_N = G u_N + g,
%
%   where u_N stands for f_N or f*_N.  Since the exact solution satisfies
%   f = G f + g, subtracting the two relations yields the identity
%
%       f - ftilde_N = G (f - u_N).
%
%   The iterated error is therefore obtained by applying G directly to the
%   collocation error e = f - u_N.  This is the numerically stable route.
%   Forming G u_N + g and then subtracting f differences two quantities of
%   order one, so that the quadrature error of the order-one integral becomes
%   an absolute floor on the result (visible as a stalled EOC once the true
%   error drops below it).  Integrating e instead makes the quadrature error
%   relative to the small quantity being measured.
%
%   Inputs
%     sol                         : collocation or modified collocation
%                                   solution returned by solve_vie
%     x                           : evaluation points
%     fRef                        : function handle for f (exact solution,
%                                   or a fine reference)
%     kappaFun, alphaFun, betaFun : data of the operator G
%     nq                          : number of quadrature nodes per piece
%                                   (default 32)
%
%   The integral is split at the mesh points of sol, across which the error
%   is only piecewise smooth, according to the following rules.
%
%     sigma_1                : Gauss--Jacobi with the frozen weight s^(-bfix),
%                              where bfix is the maximum of beta.  The residual
%                              factor s^(bfix-beta(x)) is carried in the
%                              integrand.
%     sigma_k, 2 <= k <= n-1 : Gauss--Legendre on fixed nodes.
%     [s_{n-1}, x]           : Gauss--Jacobi with weight (x-s)^(-alpha(x)).
%     n = 1                  : Gauss--Jacobi with weight
%                              (x-s)^(-alpha(x)) s^(-beta(x)).

    if nargin < 7 || isempty(nq), nq = 32; end
    if isempty(betaFun), betaFun = @(t) zeros(size(t)); end

    t = sol.t;  h = sol.h;  N = numel(h);  T = t(end);
    x = x(:);

    eStar = @(s) fRef(s) - sol.eval(s);

    bfix = max(betaFun(linspace(0,T,201).'));
    gjCache = containers.Map('KeyType','char','ValueType','any');

    [uJ0,wJ0] = getGJ(gjCache,nq,0,-bfix);
    [uL ,wL ] = getGJ(gjCache,nq,0,0);
    sFix = cell(N,1);
    sFix{1} = h(1)*(uJ0+1)/2;
    for k = 2:N
        sFix{k} = t(k) + h(k)*(uL+1)/2;
    end
    allS = vertcat(sFix{:});
    allE = eStar(allS);
    eFix = cell(N,1);  off = 0;
    for k = 1:N
        eFix{k} = allE(off+(1:nq));  off = off + nq;
    end

    P = numel(x);
    nOf = zeros(P,1);  uw = cell(P,1);  sLoc = cell(P,1);
    for q = 1:P
        xq = x(q);
        if xq < 0 || xq > T, error('eval point outside [0,T].'); end
        if xq == 0, continue; end

        n = find(t >= xq,1,'first') - 1;
        if n < 1, n = 1; end
        if n > N, n = N; end
        if xq - t(n) <= 0, n = max(n-1,1); end
        nOf(q) = n;
        ax = alphaFun(xq);  bx = betaFun(xq);
        if n == 1
            [u,w] = getGJ(gjCache,nq,-ax,-bx);
            sLoc{q} = xq*(u+1)/2;
        else
            [u,w] = getGJ(gjCache,nq,-ax,0);
            sLoc{q} = t(n) + (xq-t(n))*(u+1)/2;
        end
        uw{q} = w;
    end
    act  = find(nOf > 0);
    eLoc = cell(P,1);
    if ~isempty(act)
        allE = eStar(vertcat(sLoc{act}));
        off = 0;
        for q = act.'
            eLoc{q} = allE(off+(1:nq));  off = off + nq;
        end
    end

    vals = zeros(P,1);
    for q = act.'
        xq = x(q);  n = nOf(q);  w = uw{q};  s = sLoc{q};
        ax = alphaFun(xq);  bx = betaFun(xq);
        if n == 1
            I  = (xq/2)^(1-ax-bx) * sum(w .* kappaFun(xq,s) .* eLoc{q});
        else
            L = xq - t(n);
            I = (L/2)^(1-ax) * sum(w .* s.^(-bx) .* kappaFun(xq,s) .* eLoc{q});
            s = sFix{1};
            I = I + (h(1)/2)^(1-bfix) * ...
                sum(wJ0 .* (xq-s).^(-ax) .* s.^(bfix-bx) .* kappaFun(xq,s) .* eFix{1});
            for k = 2:n-1
                s = sFix{k};
                I = I + (h(k)/2) * ...
                    sum(wL .* (xq-s).^(-ax) .* s.^(-bx) .* kappaFun(xq,s) .* eFix{k});
            end
        end
        vals(q) = abs(I);
    end
    err = max(vals);
end

function [u,w] = getGJ(cacheMap,nq,a,b)
    key = sprintf('%.12g_%.12g',a,b);
    if isKey(cacheMap,key)
        v = cacheMap(key);  u = v{1};  w = v{2};
    else
        [u,w] = gauss_jacobi(nq,a,b);
        u = u(:);  w = w(:);
        cacheMap(key) = {u,w};
    end
end
