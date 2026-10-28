function sol = solve_vie(gFun,kappaFun,alphaFun,betaFun,T,N,m,r,c,opts)
%SOLVE_VIE  Collocation and modified collocation for a doubly weakly singular VIE.
%
%   The function solves
%
%       f(s) - int_0^s (s-t)^(-alpha(s)) t^(-beta(s)) kappa(s,t) f(t) dt = g(s)
%
%   in the space S_{m-1}^{(-1)} on the graded mesh s_n = T (n/N)^r.  Both
%   methods are marched from one subinterval to the next and both require O(N^2) kernel evaluations.
%
%   The modified solution is written as f*_N = phi_N + z_N, where
%
%       w_N = g + G phi_N,   z_N = (I - P_N) w_N,
%       phi_N(s_ni) = w_N(s_ni) + (G z_N)(s_ni).
%
%   This is the piecewise modified equation evaluated at the collocation points.
%   Once a subinterval is solved, z_N is stored at the quadrature nodes of that subinterval,
%   so that G z_N at later collocation points reduces to a quadrature sum.
%   Neither the matrix C_N nor any double integral is formed.
%
%   Inputs
%     gFun, kappaFun(s,t), alphaFun, betaFun : data of the equation
%                                              (betaFun may be empty)
%     T, N, m, r                             : interval, number of subintervals,
%                                              polynomial degree plus one,
%                                              and grading exponent
%     c                                      : collocation parameters
%                                              ([] selects Gauss points)
%     opts.method : 'colloc' (default) or 'modified'
%     opts.Q      : number of Gauss nodes per quadrature piece (default 16)
%     opts.graded : true (default).  Subintervals adjacent to the origin are split
%                   geometrically so that t^(-beta(s)) is resolved on strongly
%                   graded meshes.  If false, a single Gauss rule is used on
%                   each subinterval.
%
%   Output
%     sol.t, sol.h, sol.c, sol.Y   : mesh, steps, collocation parameters, and
%                                    values of phi_N (or f_N) at the collocation
%                                    points
%     sol.wcol                     : values of w_N at the collocation points
%     sol.eval                     : f_N (collocation) or f*_N (modified)
%     sol.eval_w                   : w_N = g + G phi_N (iterated collocation
%                                    solution when the method is collocation)
%     sol.eval_phi                 : the piecewise polynomial phi_N

    if nargin < 10 || isempty(opts), opts = struct(); end
    method = 'colloc';  if isfield(opts,'method'), method = opts.method; end
    Q      = 16;        if isfield(opts,'Q'),      Q      = opts.Q;      end
    graded = true;      if isfield(opts,'graded'), graded = opts.graded; end
    modified = strcmpi(method,'modified');
    if isempty(betaFun), betaFun = @(s) zeros(size(s)); end
    hasBeta = any(betaFun(linspace(0,T,201)) ~= 0);
    split0  = graded && hasBeta;

    if isempty(c)
        xc = gauss_jacobi(m,0,0);
        c  = (xc(:).' + 1)/2;
    end
    c = sort(c(:).');  m = numel(c);
    if any(c <= 0) || any(c > 1), error('collocation parameters must lie in (0,1].'); end

    t = T*((0:N)/N).^r;  h = diff(t);

    [xg,wg] = gauss_jacobi(Q,0,0);
    xg = (xg(:)+1)/2;  wg = wg(:)/2;
    gj = containers.Map('KeyType','char','ValueType','any');

    S = cell(N,1);  W = cell(N,1);  Lf = cell(N,1);
    first = zeros(N,1);  last = zeros(N,1);  cnt = 0;
    for k = 2:N
        [S{k},W{k}] = cell_rule(t(k),t(k+1));
        Lf{k} = lagr((S{k}-t(k))/h(k));
        first(k) = cnt+1;  cnt = cnt + numel(S{k});  last(k) = cnt;
    end
    if N >= 2, allS = vertcat(S{2:N});  allW = vertcat(W{2:N});
    else,      allS = zeros(0,1);       allW = zeros(0,1);       end
    PH = zeros(cnt,1);
    Z  = zeros(cnt,1);

    Y    = zeros(N,m);
    wcol = zeros(N,m);

    for n = 1:N
        X  = t(n) + c(:)*h(n);
        gv = col(gFun(X));
        if n == 1
            H = zeros(m,1);  Zh = zeros(m,1);
        else
            [H,Zh] = hist_at(X,n,modified);
        end
        A = zeros(m,m);  tau = cell(m,1);  wt = cell(m,1);
        for i = 1:m
            [A(i,:),tau{i},wt{i}] = cur_rule(X(i),n);
        end

        if ~modified
            a = (eye(m) - A) \ (gv + H);
        else
            TT = vertcat(tau{:});
            if n == 1, HT = zeros(size(TT)); else, HT = hist_at(TT,n,false); end
            RT = cur_rows(TT,n);
            LT = lagr((TT - t(n))/h(n));
            z0 = col(gFun(TT)) + HT - LT*(gv + H);
            Zc = RT - LT*A;
            B = zeros(m,m);  rz = zeros(m,1);  off = 0;
            for i = 1:m
                ii = off + (1:numel(tau{i}));  off = off + numel(tau{i});
                B(i,:) = wt{i}.' * Zc(ii,:);
                rz(i)  = wt{i}.' * z0(ii);
            end
            a = (eye(m) - A - B) \ (gv + H + Zh + rz);
        end

        Y(n,:)    = a.';
        wcol(n,:) = (gv + H + A*a).';
        if n >= 2 && n <= N-1
            ii = first(n):last(n);
            PH(ii) = Lf{n}*a;
            if modified
                Sn = S{n};
                wn = col(gFun(Sn)) + hist_at(Sn,n,false) + cur_rows(Sn,n)*a;
                Z(ii) = wn - Lf{n}*wcol(n,:).';
            end
        end
    end

    sol.t = t;  sol.h = h;  sol.c = c;  sol.Y = Y;  sol.wcol = wcol;
    sol.method = method;  sol.Q = Q;  sol.graded = graded;
    sol.nodes = cnt;
    sol.eval_phi = @(x) eval_parts(x,1);
    sol.eval_w   = @(x) eval_parts(x,2);
    if modified
        sol.eval = @(x) eval_parts(x,3);
    else
        sol.eval = sol.eval_phi;
    end

    function [s,w] = cell_rule(a,b)
        bp = a;
        if split0
            p = a;
            while 2*p < b, p = 2*p;  bp(end+1) = p; end %#ok<AGROW>
        end
        bp(end+1) = b;
        np = numel(bp)-1;
        s = zeros(Q*np,1);  w = s;
        for j = 1:np
            L = bp(j+1) - bp(j);
            s((j-1)*Q+(1:Q)) = bp(j) + L*xg;
            w((j-1)*Q+(1:Q)) = L*wg;
        end
    end

    function [H,Zh] = hist_at(X,n,wantZ)
        X  = col(X);  P = numel(X);
        ax = expand(alphaFun(X),P);  bx = expand(betaFun(X),P);
        d  = X - t(n);
        H  = zeros(P,1);  Zh = zeros(P,1);
        y1 = Y(1,:).';

        for q = 1:P
            [u,w] = rule(0,-bx(q));
            near1 = (n == 2) && d(q) < h(1);
            if near1, b1 = h(1)/2; else, b1 = h(1); end
            s1 = b1*(u+1)/2;
            v  = (b1/2)^(1-bx(q)) * w .* (X(q)-s1).^(-ax(q)) .* kap(X(q),s1);
            H(q) = v.' * (lagr(s1/h(1))*y1);
            if near1
                H(q) = H(q) + adjacent(X(q),ax(q),bx(q),h(1)/2,h(1),1);
            end
            if wantZ
                if near1
                    s1 = h(1)*(u+1)/2;
                    v  = (h(1)/2)^(1-bx(q)) * w .* (X(q)-s1).^(-ax(q)) .* kap(X(q),s1);
                end
                Zh(q) = v.' * z_cell1(s1);
            end
        end

        if n >= 3
            if n >= 4
                jj = 1:last(n-2);
                K  = kern(X,ax,bx,allS(jj));
                H  = H + K*(allW(jj).*PH(jj));
                if wantZ, Zh = Zh + K*(allW(jj).*Z(jj)); end
            end
            jj  = first(n-1):last(n-1);
            K   = kern(X,ax,bx,allS(jj));
            far = d >= h(n-1);
            H(far) = H(far) + K(far,:)*(allW(jj).*PH(jj));
            if wantZ, Zh = Zh + K*(allW(jj).*Z(jj)); end
            if any(~far)
                H(~far) = H(~far) + adjacent(X(~far),ax(~far),bx(~far),t(n-1),t(n),n-1);
            end
        end
    end

    function v = adjacent(Xv,axv,bxv,a,b,k)
        Xv = col(Xv);  v = zeros(size(Xv));
        yk = Y(k,:).';  t0 = t(k);  Hk = h(k);
        bp = a;
        if split0 && a > 0
            p = a;
            while 2*p < b, p = 2*p;  bp(end+1) = p; end %#ok<AGROW>
        end
        bp(end+1) = b;
        for j = 1:numel(bp)-1
            p = bp(j);  q = bp(j+1);  L = q - p;
            dd = Xv - q;  lo = dd;  act = true(size(Xv));
            while any(act)
                ia  = find(act);
                hi = min(2*lo(ia), dd(ia) + L);
                len = hi - lo(ia);
                R  = lo(ia) + len.*xg.';
                Sx = Xv(ia) - R;
                K  = R.^(-axv(ia)) .* Sx.^(-bxv(ia)) .* kappaFun(Xv(ia).*ones(1,Q),Sx);
                ph = reshape(lagr((Sx(:)-t0)/Hk)*yk, size(Sx));
                v(ia) = v(ia) + sum(len.*wg.'.*K.*ph, 2);
                lo(ia) = hi;
                act(ia) = hi < dd(ia) + L*(1 - 1e-14);
            end
        end
    end

    function z = z_cell1(s)
        s = col(s);
        wv = col(gFun(s)) + cur_rows(s,1)*Y(1,:).';
        z  = wv - lagr(s/h(1))*wcol(1,:).';
    end

    function [row,tau,wt] = cur_rule(x,n)
        a0 = alphaFun(x);  b0 = betaFun(x);
        if n == 1
            [u,w] = rule(-a0,-b0);
            tau = x*(u+1)/2;
            wt  = (x/2)^(1-a0-b0) * w .* kap(x,tau);
        else
            p = t(n);  tau = zeros(0,1);  wt = zeros(0,1);
            if split0
                while p < x/2
                    q  = min(2*p, x/2);
                    tt = p + (q-p)*xg;
                    tau = [tau; tt];                                       %#ok<AGROW>
                    wt  = [wt; (q-p)*wg.*(x-tt).^(-a0).*tt.^(-b0).*kap(x,tt)]; %#ok<AGROW>
                    p = q;
                end
            end
            [u,w] = rule(-a0,0);
            L  = x - p;  tt = p + L*(u+1)/2;
            tau = [tau; tt];
            wt  = [wt; (L/2)^(1-a0) * w .* tt.^(-b0) .* kap(x,tt)];
        end
        row = wt.' * lagr((tau - t(n))/h(n));
    end

    function R = cur_rows(Xv,n)
        Xv = col(Xv);  R = zeros(numel(Xv),m);
        for q = 1:numel(Xv), R(q,:) = cur_rule(Xv(q),n); end
    end

    function y = eval_parts(x,what)
        sz = size(x);  x = col(x);  P = numel(x);  y = zeros(P,1);
        if any(x < 0) || any(x > T*(1+1e-14)), error('eval point outside [0,T].'); end
        nx = ceil(N*(x/T).^(1/r));  nx = min(max(nx,1),N);
        i1 = nx > 1 & x <= t(nx).';     nx(i1) = nx(i1) - 1;
        i2 = nx < N & x >  t(nx+1).';   nx(i2) = nx(i2) + 1;
        u  = (x - t(nx).')./h(nx).';
        Lx = lagr(u);
        phi = sum(Lx.*Y(nx,:),2);
        if what == 1, y = reshape(phi,sz); return; end
        wv = col(gFun(x));
        for nc = unique(nx).'
            sel = find(nx == nc & x > 0);
            if isempty(sel), continue; end
            Xs = x(sel);
            if nc >= 2, Hs = hist_at(Xs,nc,false); else, Hs = 0; end
            wv(sel) = wv(sel) + Hs + cur_rows(Xs,nc)*Y(nc,:).';
        end
        if what == 2
            y = wv;
        else
            y = phi + wv - sum(Lx.*wcol(nx,:),2);
        end
        y = reshape(y,sz);
    end

    function K = kern(X,ax,bx,Sv)
        Sm = Sv(:).';
        K  = exp(-ax.*log(X - Sm) - bx.*log(Sm)) .* kappaFun(X.*ones(size(Sm)), ones(size(X)).*Sm);
    end

    function k = kap(x,tau)
        k = kappaFun(x*ones(size(tau)),tau);
    end

    function L = lagr(u)
        u = u(:);  L = ones(numel(u),m);
        for j = 1:m
            for kk = [1:j-1, j+1:m]
                L(:,j) = L(:,j) .* (u - c(kk)) / (c(j) - c(kk));
            end
        end
    end

    function [u,w] = rule(a,b)
        key = sprintf('%.17g|%.17g',a,b);
        if isKey(gj,key)
            v = gj(key);  u = v{1};  w = v{2};
        else
            [u,w] = gauss_jacobi(Q,a,b);
            u = u(:);  w = w(:);
            gj(key) = {u,w};
        end
    end
end

function v = col(v)
    v = v(:);
end

function v = expand(v,P)
    v = v(:);
    if numel(v) == 1, v = v*ones(P,1); end
end
