function info = iqc_analysis(iqc,G,K,CH,options)
if nargin == 3
    options = struct();
end
if ~isfield(options,"rho_start")
    options.rho_start = 0.99;
end
if ~isfield(options,"margin_ana") % margin for strict inequalities, i.e.,
    options.margin_ana = 5e-8;    % LMI > 0   <==>   LMI >= margin_pos
end
for j = 1:length(CH)
    optsj=options;
    CH(j).gam = -1;
    if CH(j).gain ~= "p2p"
        optsj.rho_start=1;
    end
    info(j) = ana_step(iqc,G,K,CH(j),[],optsj);
end

