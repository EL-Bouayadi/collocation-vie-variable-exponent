function plot_cost(tag,rC,rM,outName)
%PLOT_COST  Work-precision figures of the manuscript.
%
%   plot_cost
%       Default figure (m = 2).
%   plot_cost('fig1',2.5,1.25,'cost_ex1_m1')
%       Figure for m = 1.
%
%   The function reads results_cost_<tag>.mat, as produced by
%   cost_experiment({tag}), and writes ../figures/<outName>.pdf.  The plot
%   shows -log10(E_N), with E_N measured on S_N including the right limits,
%   against the CPU time (solution plus evaluation on S_N) on a logarithmic
%   scale.

    if nargin < 1
        tag = 'fig';  rC = 5;  rM = 2.5;  outName = 'cost_ex1';
    end
    thisDir = fileparts(mfilename('fullpath'));
    S = load(fullfile(thisDir,['results_cost_' tag '.mat']));  R = S.R;
    outDir = fullfile(thisDir,'..','figures');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    sel = { {'colloc',rC}, {'modified',rM}, {'modified',rC} };
    col = {[0x2a 0x78 0xd6], [0xeb 0x68 0x34], [0x1b 0xaf 0x7a]};
    mk  = {'o','s','^'};
    ls  = {'-','--','-.'};
    lab = {sprintf('$f_N$, $r_{\\mathrm C}=%g$',rC), sprintf('$f_N^*$, $r_{\\mathrm M}=%g$',rM), ...
           sprintf('$f_N^*$, $r_{\\mathrm C}=%g$',rC)};
    ink = [0.32 0.32 0.31];

    fig = figure('Units','centimeters','Position',[2 2 6.2 5.2],'Color','w');
    ax = axes(fig);  hold(ax,'on');
    h = gobjects(1,3);  dall = [];
    for k = 1:3
        i = find(strcmp({R.method},sel{k}{1}) & [R.r] == sel{k}{2}, 1);
        t = R(i).tS + R(i).tE;  d = -log10(R(i).Ep);  dall = [dall d]; %#ok<AGROW>
        c = double(col{k})/255;
        h(k) = plot(ax,t,d,'LineStyle',ls{k},'Marker',mk{k},'Color',c, ...
                    'LineWidth',1.2,'MarkerSize',4.5,'MarkerFaceColor',c,'MarkerEdgeColor','w');
    end
    set(ax,'XScale','log','YScale','linear','FontSize',7,'Box','off', ...
        'XColor',ink,'YColor',ink,'TickLabelInterpreter','latex', ...
        'XGrid','on','YGrid','on','XMinorGrid','off','YMinorGrid','off', ...
        'GridColor',[0 0 0],'GridAlpha',0.10,'LineWidth',0.5,'TickDir','out');
    y0 = floor(min(dall));  y1 = ceil(max(dall));
    ylim(ax,[y0 y1]);  yticks(ax,y0:y1);
    xl = xlim(ax);  xticks(ax,10.^(floor(log10(xl(1))):ceil(log10(xl(2)))));
    xlabel(ax,'CPU time (s)','Interpreter','latex','FontSize',8);
    ylabel(ax,'$-\log_{10}E_N$','Interpreter','latex','FontSize',8);
    legend(ax,h,lab,'Interpreter','latex','FontSize',6.5,'Location','northwest', ...
           'Box','on','EdgeColor',ink,'LineWidth',0.5);

    out = fullfile(outDir,[outName '.pdf']);
    exportgraphics(fig,out,'ContentType','vector');
    exportgraphics(fig,fullfile(outDir,[outName '.png']),'Resolution',200);
    fprintf('wrote %s\n', out);
end
