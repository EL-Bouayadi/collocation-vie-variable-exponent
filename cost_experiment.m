function cost_experiment(which)
%COST_EXPERIMENT  CPU time and accuracy for the cost figures of the paper.
%
%   For each example and each method the function records the solve time, the
%   time to evaluate the approximation on the sample set S_N of the manuscript
%   (21 points per subinterval), and the errors E_N (on S_N) and E_N^+ (on S_N together
%   with the right limits at the left end of every subinterval, where the error of a
%   discontinuous approximation is often largest).  All solves use solve_vie,
%   so that collocation and modified collocation share the same quadrature.
%
%   Methods
%     colloc    f_N         on r_C (and on r_M)
%     iter      g + G f_N   on r_M (iterated collocation, any c)
%     modified  f*_N        on r_M (and on r_C)
%
%   Example A is Example 1 of the manuscript (gamma_0 = 0.6), with m = 2 and
%   c = (1/2,1).  Example B has gamma_0 = 0.9, alpha(s) = 0.5-0.2s,
%   beta(s) = 0.4-0.2s and f(s) = 1 + s^0.1, with m = 1 (c = 1) and m = 2
%   (c = (1/2,1)).  Example B also runs collocation on r_C with a single Gauss
%   rule per subinterval (graded = false), in order to exhibit the failure of that
%   naive quadrature.
%
%   The recommended call is single-threaded:
%     matlab -singleCompThread -batch "cost_experiment"

    if nargin < 1, which = {'A','B1','B2'}; end
    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);
    tag = '';
    if ~isequal(which,{'A','B1','B2'}), tag = ['_' strjoin(which,'')]; end
    outFile = fullfile(thisDir,['log_cost' tag '.txt']);
    if exist(outFile,'file'), delete(outFile); end
    diary(outFile);  cl = onCleanup(@() diary('off')); %#ok<NASGU>

    fprintf('COST_EXPERIMENT  %s   MATLAB %s   threads %d\n', datestr(now), version, maxNumCompThreads);
    R = struct('ex',{},'m',{},'method',{},'r',{},'graded',{},'N',{},'E',{},'Ep',{},'tS',{},'tE',{},'qc',{});

    for w = which
        switch w{1}
            case 'A'
                d = data_ex(0.45,-0.25,0.15,-0.10);
                m = 2;  c = [1/2 1];  rC = 5;  rM = 2.5;
                R = [R, run_block('A',d,m,c,'colloc',  rC,true,[16 32 64 128 256 512 1024 2048])]; %#ok<AGROW>
                R = [R, run_block('A',d,m,c,'colloc',  rM,true,[16 32 64 128 256 512 1024 2048])]; %#ok<AGROW>
                R = [R, run_block('A',d,m,c,'iter',    rM,true,[16 32 64 128 256 512 1024])];      %#ok<AGROW>
                R = [R, run_block('A',d,m,c,'modified',rM,true,[16 32 64 128 256 512 1024])];      %#ok<AGROW>
                R = [R, run_block('A',d,m,c,'modified',rC,true,[16 32 64 128 256 512])];           %#ok<AGROW>
            case {'B1','B2'}
                d = data_ex(0.5,-0.2,0.4,-0.2);
                if strcmp(w{1},'B1'), m = 1; c = 1; else, m = 2; c = [1/2 1]; end
                rC = m/0.1;  rM = rC/2;
                NL = [16 32 64 128 256 512 1024];
                R = [R, run_block(w{1},d,m,c,'colloc',  rC,true, [NL 2048])]; %#ok<AGROW>
                R = [R, run_block(w{1},d,m,c,'colloc',  rC,false,NL)];        %#ok<AGROW>
                R = [R, run_block(w{1},d,m,c,'iter',    rM,true, NL)];        %#ok<AGROW>
                R = [R, run_block(w{1},d,m,c,'modified',rM,true, NL)];        %#ok<AGROW>
                R = [R, run_block(w{1},d,m,c,'modified',rM,false,NL)];        %#ok<AGROW>
            case {'B1C','B2C'}
                d = data_ex(0.5,-0.2,0.4,-0.2);
                if strcmp(w{1},'B1C'), m = 1; c = 1; else, m = 2; c = [1/2 1]; end
                R = [R, run_block(w{1},d,m,c,'modified',m/0.1,true,[16 32 64 128 256 512])]; %#ok<AGROW>
            case 'fig'
                d = data_ex(0.45,-0.25,0.15,-0.10);
                R = [R, run_interleaved(d,2,[1/2 1], ...
                        {{'colloc',5,2.^(4:12)},{'modified',2.5,2.^(4:10)},{'modified',5,2.^(4:10)}},3)]; %#ok<AGROW>
            case 'fig1'
                d = data_ex(0.45,-0.25,0.15,-0.10);
                R = [R, run_interleaved(d,1,1, ...
                        {{'colloc',2.5,2.^(4:12)},{'modified',1.25,2.^(4:10)},{'modified',2.5,2.^(4:10)}},3)]; %#ok<AGROW>
        end
        save(fullfile(thisDir,['results_cost' tag '.mat']),'R');
    end
    fprintf('done %s\n', datestr(now));
end

function d = data_ex(a0,a1,b0,b1)
    d.T = 1;  d.kap = -1;
    d.alphaFun = @(s) a0 + a1*s;
    d.betaFun  = @(s) b0 + b1*s;
    d.kappaFun = @(s,t) d.kap*ones(size(t));
    d.p = 1 - a0 - b0;
    d.f = @(s) 1 + s.^d.p;
    d.g = @(s) mrhs(s,d.alphaFun,d.betaFun,d.kap,d.p);
    d.label = sprintf('alpha=%g%+gs beta=%g%+gs',a0,a1,b0,b1);
end

function g = mrhs(s,alphaFun,betaFun,kap,p)
    s = s(:);  a = alphaFun(s);  b = betaFun(s);
    Gf = beta(1-a,1-b).*s.^(1-a-b) + beta(1-a,1-b+p).*s.^(1-a-b+p);
    Gf(s == 0) = 0;
    g = (1 + s.^p) - kap*Gf;
end

function res = run_interleaved(d,m,c,cfg,nsweep)
    K = numel(cfg);
    for k = 1:K
        res(k) = struct('ex','fig','m',m,'method',cfg{k}{1},'r',cfg{k}{2}, ...
                        'graded',true,'N',cfg{k}{3},'E',[],'Ep',[], ...
                        'tS',inf(size(cfg{k}{3})),'tE',inf(size(cfg{k}{3})),'qc',NaN); %#ok<AGROW>
    end
    lists = cellfun(@(x) x{3}, cfg, 'UniformOutput', false);
    allN = unique([lists{:}]);
    for sw = 1:nsweep
        fprintf('\n=== sweep %d of %d  (%s) ===\n', sw, nsweep, datestr(now));
        for N = allN
            for k = 1:K
                i = find(res(k).N == N, 1);
                if isempty(i), continue; end
                opts = struct('method',res(k).method,'graded',true);
                tic;
                sol = solve_vie(d.g,d.kappaFun,d.alphaFun,d.betaFun,d.T,N,m,res(k).r,c,opts);
                ts = toc;
                [xs,isp] = sample_points(sol.t);
                tic;  v = sol.eval(xs);  te = toc;
                res(k).tS(i) = min(res(k).tS(i), ts);
                res(k).tE(i) = min(res(k).tE(i), te);
                if sw == 1
                    err = abs(v - d.f(xs));
                    res(k).E(i) = max(err(~isp));  res(k).Ep(i) = max(err);
                end
                fprintf('  %-8s r=%-4g N=%5d  E+ %.4e  solve %8.3f  eval %7.3f\n', ...
                        res(k).method, res(k).r, N, res(k).Ep(i), ts, te);
            end
        end
    end
    fprintf('\n--- minima over %d sweeps ---\n', nsweep);
    for k = 1:K
        fprintf('%s r=%g\n%6s %12s %10s %10s %10s\n', res(k).method, res(k).r, ...
                'N','E_N^+','solve[s]','eval[s]','total[s]');
        fprintf('%6d %12.4e %10.3f %10.3f %10.3f\n', ...
                [res(k).N; res(k).Ep; res(k).tS; res(k).tE; res(k).tS+res(k).tE]);
    end
end

function res = run_block(ex,d,m,c,method,r,graded,NL,nrep)
    if nargin < 9, nrep = 1; end
    fprintf('\n--- %s  %s  m=%d c=%s  %s  r=%g  graded=%d  (min of %d timings) ---\n', ...
            ex, d.label, m, mat2str(c,3), method, r, graded, nrep);
    fprintf('%6s %12s %6s %12s %6s %10s %10s\n','N','E_N','EOC','E_N^+','EOC','solve[s]','eval[s]');
    if strcmp(method,'iter'), solver = 'colloc'; else, solver = method; end
    E = nan(size(NL));  Ep = E;  tS = E;  tE = E;  qc = NaN;
    for i = 1:numel(NL)
        N = NL(i);
        opts = struct('method',solver,'graded',graded);
        tS(i) = Inf;  tE(i) = Inf;
        for rep = 1:nrep
            tic;
            sol = solve_vie(d.g,d.kappaFun,d.alphaFun,d.betaFun,d.T,N,m,r,c,opts);
            tS(i) = min(tS(i), toc);
            [xs,isp] = sample_points(sol.t);
            tic;
            if strcmp(method,'iter'), v = sol.eval_w(xs); else, v = sol.eval(xs); end
            tE(i) = min(tE(i), toc);
        end
        err = abs(v - d.f(xs));
        E(i) = max(err(~isp));  Ep(i) = max(err);
        if i > 1
            eoc = log2(E(i-1)/E(i));  eocp = log2(Ep(i-1)/Ep(i));
        else
            eoc = NaN;  eocp = NaN;
        end
        fprintf('%6d %12.4e %6.2f %12.4e %6.2f %10.3f %10.3f\n', N, E(i), eoc, Ep(i), eocp, tS(i), tE(i));
        if i == numel(NL) || tS(i) + tE(i) > 900
            o2 = opts;  o2.Q = 32;
            s2 = solve_vie(d.g,d.kappaFun,d.alphaFun,d.betaFun,d.T,N,m,r,c,o2);
            if strcmp(method,'iter'), v2 = s2.eval_w(xs); else, v2 = s2.eval(xs); end
            qc = abs(max(abs(v2 - d.f(xs))) - Ep(i)) / Ep(i);
            fprintf('   quadrature check (Q 16 -> 32) at N=%d: rel. change of E_N %.1e\n', N, qc);
            if i < numel(NL)
                fprintf('   (stopped: time limit)\n');
                E(i+1:end) = [];  Ep(i+1:end) = [];  tS(i+1:end) = [];  tE(i+1:end) = [];  NL(i+1:end) = [];
            end
            break
        end
    end
    res = struct('ex',ex,'m',m,'method',method,'r',r,'graded',graded,'N',NL, ...
                 'E',E,'Ep',Ep,'tS',tS,'tE',tE,'qc',qc);
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
