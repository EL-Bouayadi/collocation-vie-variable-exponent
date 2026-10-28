function R = numerical_examples(ex,blocks,NL,NmaxIM)
%NUMERICAL_EXAMPLES  Global errors for Examples 1--4 of the paper.
%
%   For each pair (m,r) the function prints the global errors and the
%   estimated order of convergence on every row.
%
%   The seven columns, in the same order as the manuscript, with an extra
%   column for f*_N at Gauss points, are as follows.
%     1  f_N     with collocation parameters c for which J_0 is nonzero
%     2  ft_N    with Gauss parameters (J_0 = 0)
%     3  ft_N    with collocation parameters c for which J_0 is nonzero
%     4  f*_N    with collocation parameters c for which J_0 is nonzero
%     5  f*_N    with Gauss parameters
%     6  ft*_N   with Gauss parameters
%     7  ft*_N   with collocation parameters c for which J_0 is nonzero
%
%   The parameters c are 1, (1/2,1) and (1/3,2/3,1) for m = 1, 2 and 3.
%   Iterated errors are computed as G(f-u_N) by iterated_error.  The remaining
%   errors are maxima over the sample set S_N, which includes the right limits
%   s_{n-1} + 1e-10 h_n.
%
%   Usage
%     numerical_examples(1)
%         Example 1 with the default blocks.
%     numerical_examples(3,[2 1; 2 2.5],[8 16 32])
%         Selected blocks, given as rows [m r].
%     numerical_examples(4)
%         Example with beta_* > gamma_0.
%
%   Example 3 uses the reference solution of the paper (collocation with
%   m = 4, r = 10 and N = 1024).  The floor ||u_512 - u_1024|| is printed.

    if nargin < 2, blocks = []; end
    if nargin < 3 || isempty(NL),     NL     = [16 32 64 128 256 512]; end
    if nargin < 4 || isempty(NmaxIM), NmaxIM = 256; end

    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    outFile = fullfile(thisDir,sprintf('log_example_%d.txt',ex));
    if exist(outFile,'file'), delete(outFile); end
    diary(outFile);  cl = onCleanup(@() diary('off')); %#ok<NASGU>

    d = example_data(ex);
    if isempty(blocks), blocks = d.blocks; end
    fprintf('NUMERICAL_EXAMPLES  Example %d  %s  MATLAB %s\n', ex, datestr(now), version);
    fprintf('%s\n', d.label);

    if d.needRef
        tic;
        sRef  = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,1024,4,10,[],struct('method','colloc'));
        sHalf = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,512,4,10,[],struct('method','colloc'));
        xF = sample_points(sRef.t);
        fprintf('reference m=4 r=10 N=1024: floor ||u_512-u_1024|| = %.3e  (%.0f s)\n', ...
                max(abs(sHalf.eval(xF) - sRef.eval(xF))), toc);
        d.f = @(s) reshape(sRef.eval(s(:)),size(s));
    end

    names = {'f_N c','ft_N G','ft_N c','f*_N c','f*_N G','ft*_N G','ft*_N c'};
    R = struct('m',{},'r',{},'N',{},'E',{});
    tex = {};
    for b = 1:size(blocks,1)
        m = blocks(b,1);  r = blocks(b,2);
        c = cJ0(m);
        th = theory(d,m,r);
        fprintf('\n--- m = %d   r = %.4g   c = %s\n', m, r, mat2str(c,4));
        fprintf('theory: muC %.2f  muI(G) %.2f  muI0(c) %.2f  muM %.2f  muIM %.2f\n', th);
        fprintf('%6s', 'N');  fprintf(' %11s %5s', names{:});  fprintf(' %8s\n','sec');
        E = nan(numel(NL),7);
        for i = 1:numel(NL)
            N = NL(i);  tic;
            sC  = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,N,m,r,c, struct('method','colloc'));
            sG  = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,N,m,r,[],struct('method','colloc'));
            sM  = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,N,m,r,c, struct('method','modified'));
            sMG = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,N,m,r,[],struct('method','modified'));
            [xs,isp] = sample_points(sC.t);
            fx = d.f(xs);  xc = xs(~isp);
            E(i,1) = max(abs(fx - sC.eval(xs)));
            E(i,2) = iterated_error(sG,xc,d.f,d.kappa,d.alpha,d.beta,32);
            E(i,3) = iterated_error(sC,xc,d.f,d.kappa,d.alpha,d.beta,32);
            E(i,4) = max(abs(fx - sM.eval(xs)));
            E(i,5) = max(abs(fx - sMG.eval(xs)));
            if N <= NmaxIM
                E(i,6) = iterated_error(sMG,xc,d.f,d.kappa,d.alpha,d.beta,32);
                E(i,7) = iterated_error(sM, xc,d.f,d.kappa,d.alpha,d.beta,32);
            end
            eoc = nan(1,7);
            if i > 1, eoc = log2(E(i-1,:)./E(i,:)); end
            fprintf('%6d', N);
            for k = 1:7, fprintf(' %11.3e %5.2f', E(i,k), eoc(k)); end
            fprintf(' %8.1f\n', toc);
            R(end+1) = struct('m',m,'r',r,'N',N,'E',E(i,:)); %#ok<AGROW>
            save(fullfile(thisDir,sprintf('results_example_%d.mat',ex)),'R');
        end
        N = NL(end);  o2 = struct('method','colloc','Q',32);
        sC2 = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,N,m,r,c,o2);
        [xs,isp] = sample_points(sC2.t);
        q1 = abs(max(abs(d.f(xs) - sC2.eval(xs))) - E(end,1))/E(end,1);
        q3 = abs(iterated_error(sC,xs(~isp),d.f,d.kappa,d.alpha,d.beta,64) - E(end,3))/E(end,3);
        o2.method = 'modified';
        sM2 = solve_vie(d.g,d.kappa,d.alpha,d.beta,1,N,m,r,c,o2);
        q4 = abs(max(abs(d.f(xs) - sM2.eval(xs))) - E(end,4))/E(end,4);
        fprintf('   quadrature check at N=%d: f_N %.1e  ft_N %.1e  f*_N %.1e (relative change)\n', N, q1, q3, q4);
        tex = [tex, latex_block(m,r,d,NL,E,th)]; %#ok<AGROW>
    end
    texFile = fullfile(thisDir,sprintf('tables_example_%d.tex',ex));
    fid = fopen(texFile,'w');  fprintf(fid,'%s\n',tex{:});  fclose(fid);
    fprintf('\nLaTeX rows written to %s\ndone %s\n', texFile, datestr(now));
end

function d = example_data(ex)
    d.needRef = false;
    switch ex
        case 1
            d.alpha = @(s) 0.45 - 0.25*s;  d.beta = @(s) 0.15 - 0.10*s;
            d.kappa = @(s,t) -ones(size(t));
            d.f = @(s) 1 + s.^0.4;
            d.g = @(s) rhs_power(s,d.alpha,d.beta,0.4);
            d.g0 = 0.6;  d.as = 0.45;  d.bs = 0.15;
            d.blocks = [1 1; 1 1.25; 1 2.5; 2 1; 2 2.5; 2 5];
            d.label = 'Example 1: alpha=0.45-0.25s beta=0.15-0.10s f=1+s^0.4';
        case 2
            d.alpha = @(s) 0.6*s.^3.*(4-3*s);  d.beta = d.alpha;
            d.kappa = @(s,t) -ones(size(t));
            d.f = @(s) 1 + s + s.^2 + s.^3;
            d.g = @(s) rhs_poly(s,d.alpha,d.beta);
            d.g0 = 0;  d.as = 0.6;  d.bs = 0.6;
            d.blocks = [1 1; 2 1; 3 1];
            d.label = 'Example 2: alpha=beta=0.6 s^3(4-3s) f=1+s+s^2+s^3 (smooth case)';
        case 3
            a0 = 0.6;  a1 = 0.4;
            d.alpha = @(s) a1 + (a0-a1)*(1 - s - sin(2*pi*(1-s))/(2*pi));
            d.beta  = [];
            d.kappa = @(s,t) -ones(size(t))./gamma(1-d.alpha(s));
            d.g = @(s) -ones(size(s));
            d.f = [];  d.needRef = true;
            d.g0 = 0.6;  d.as = 0.6;  d.bs = 0;
            d.blocks = [1 1; 1 1.25; 1 2.5; 2 1; 2 2.5; 2 5];
            d.label = 'Example 3: Liang-Stynes test problem, alpha(0)=0.6 alpha(1)=0.4, beta=0';
        case 4
            d.alpha = @(s) 0.3*ones(size(s));  d.beta = @(s) 0.1 + 0.6*s;
            d.kappa = @(s,t) -ones(size(t));
            d.f = @(s) 1 + s.^0.6;
            d.g = @(s) rhs_power(s,d.alpha,d.beta,0.6);
            d.g0 = 0.4;  d.as = 0.3;  d.bs = 0.7;
            d.blocks = [2 1; 2 5/3; 2 10/3];
            d.label = 'Example 4: alpha=0.3 beta=0.1+0.6s (beta_*=0.7 > gamma_0=0.4) f=1+s^0.6';
    end
end

function th = theory(d,m,r)
    g0 = d.g0;  as = d.as;  bs = d.bs;
    if g0 == 0 && ex2like(d)
        muC = m;  muM = m+1-as;  muI = min(m+1-as, r*(m+1-bs));  muI0 = m;  muIM = m+1-as;
    else
        muC = min(r*(1-g0),m);  muM = min(2*r*(1-g0),m+1-as);
        muI = min(muM, r*(2-g0-bs));  muI0 = min([2*r*(1-g0), r*(2-g0-bs), m]);
        nu = min(r*(1-g0),1-as);
        muIM = min([muM+nu, muC+m, r*(3-2*g0-bs), m+1-as]);
    end
    th = [muC muI muI0 muM muIM];
end

function t = ex2like(d)
    t = abs(d.alpha(0)) < 1e-14;
end

function c = cJ0(m)
    switch m
        case 1, c = 1;
        case 2, c = [1/2 1];
        case 3, c = [1/3 2/3 1];
    end
end

function L = latex_block(m,r,d,NL,E,th)
    if abs(r-1) < 1e-12, rl = '$r=1$';
    else, rl = sprintf('$r=%s$', strtrim(sprintf('%.4g',r))); end
    L = {sprintf('%% m=%d r=%g  (%s)', m, r, d.label)};
    for i = 1:numel(NL)
        row = sprintf(' & $%d$', NL(i));
        for k = 1:7
            if isnan(E(i,k)), row = [row ' & --']; continue; end %#ok<AGROW>
            s = sprintf('%.2e', E(i,k));  s = regexprep(s,'e([+-])0*(\d)','e$1$2');
            if i > 1
                row = [row sprintf(' & %s\\,{\\scriptsize(%.2f)}', s, log2(E(i-1,k)/E(i,k)))]; %#ok<AGROW>
            else
                row = [row ' & ' s]; %#ok<AGROW>
            end
        end
        if i == 1
            row = [sprintf('\\multirow{%d}{*}{\\shortstack{$m=%d$\\\\ %s}}', numel(NL)+1, m, rl) row]; %#ok<AGROW>
        end
        L{end+1} = [row ' \\']; %#ok<AGROW>
    end
    L{end+1} = sprintf(' & Theory & $%.3g$ & $%.3g$ & $%.3g$ & $%.3g$ & $%.3g$ & $%.3g$ & $%.3g$ \\\\', ...
                       th(1), th(2), th(3), th(4), th(4), th(5), th(5));
    L{end+1} = '\midrule';
end

function g = rhs_power(s,alphaFun,betaFun,p)
    s = s(:);  a = alphaFun(s);  b = betaFun(s);
    Gf = beta(1-a,1-b).*s.^(1-a-b) + beta(1-a,1-b+p).*s.^(1-a-b+p);
    Gf(s == 0) = 0;
    g = (1 + s.^p) + Gf;
end

function g = rhs_poly(s,alphaFun,betaFun)
    s = s(:);  a = alphaFun(s);  b = betaFun(s);
    Gf = zeros(size(s));
    for k = 0:3, Gf = Gf + beta(1-a,1-b+k).*s.^(1-a-b+k); end
    Gf(s == 0) = 0;
    g = (1 + s + s.^2 + s.^3) + Gf;
end

function [xs,isp] = sample_points(t)
    N = numel(t)-1;  h = diff(t);
    u = (0:20)/20;
    X = t(1:N).' + h(:)*u;
    xs = X(:);  xs(xs == 0) = [];
    xp = t(1:N).' + 1e-10*h(:);
    isp = [false(numel(xs),1); true(N,1)];
    xs = [xs; xp];
end
