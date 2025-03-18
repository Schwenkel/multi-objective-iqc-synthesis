function [K,ana_ws,info] = iqc_synthesis(G,CH,iqc,options)
% CH is struct array. 
% CH(i).w contains performance input G.InputGroup.w(CH(i).w)
% CH(i).z contains performance output G.OutputGroup.z(CH(i).z)
% CH(i).gain contains performance measure "p2p", "e2p", or "e2e".
% CH(i).gam constains -1 if gam shall be optimized or a fixed positive
%           value if gain shall be less than gam. 
% set default options if not specified
if nargin == 3
    options = struct();
end
if ~isfield(options,"iterations") % number of analysis/synthesis steps 
    options.iterations = 10;      % (a/s-step)
end
if ~isfield(options,"tau_iter")   % number of iterations over tau per
    options.tau_iter = 5;         % a/s-step if tau=1 is not feasible
end
if ~isfield(options,"ind_ana_each_step")  % if true an indiviual analysis
    options.ind_ana_each_step = false;    % of each performance channel is
end                                       % performed after each a-step
if ~isfield(options,"rho_start")  % start rho optimization in each a-step 
    options.rho_start = 1;        % for p2p with rho_start, else rho = 1
    for j = 1:length(CH)
        if CH(j).gain == "p2p"
            options.rho_start = 0.99;
             break;
        end
    end
end
% warm start with nominal controller:
tic
[K, info_nms] = nom_synthesis(G,CH);
% Start iteration
info = struct;
info.K = K;
info.tau_s = 0;
info.gam_s = info_nms.gam;
ana_ws = []; % start without warm start
tau_iter = options.tau_iter;
tau_min = 0;
for i=1:options.iterations
    % analysis step
    tau_max = 1; tau = 1; solfound = false;  % initialization
    for j = 1:tau_iter
        Gt = scale_Delta_by_tau(G, tau);
        [ana_ws_try,syn_init_try] = ana_step(iqc,Gt,K,CH,ana_ws,options);
        solved = (ana_ws_try.info == "solution found");
        if solved
            solfound = true;
            syn_init = syn_init_try;
        end
        [tau,tau_min,tau_max] = update_tau(tau,tau_min,tau_max,j,solved);
        if solved && tau == 1
            break;
        end
    end
    % return if no solution found
    if ~solfound
        warning("Iteration stopped at after "+num2str(i-1)+...
                               " steps due to numerical problems.")
        if i>1
            return
        else
            error("No solution found.")
        end
    end
    info(i).gam = sum(syn_init.gam);
    info(i).tau = tau_min;
    info(i).rho = syn_init.rho;
    info(i).time = toc;
    if options.ind_ana_each_step
        info(i).ind_ana = iqc_analysis(iqc,G,info(i).K,CH,options);
    end
    display_info(info)
    tic
    % synthesis step
    tau_max = 1; tau = 1; solfound = false;  % initialization
    for j = 1:tau_iter
        Gt = scale_Delta_by_tau(G, tau);
        [K_try, ana_ws_try] = syn_step(Gt,CH,syn_init,options);
        solved = (ana_ws_try.info == "solution found");
        if solved
            solfound = true;
            K = K_try;
            ana_ws = ana_ws_try;
        end
        [tau,tau_min,tau_max] = update_tau(tau,tau_min,tau_max,j,solved);
        if solfound && tau == 1
            break;
        end
    end
    info(i+1).K = K;
    info(i+1).tau_s = tau_min;
    info(i+1).gam_s = sum(ana_ws.gam);
    % return if no solution found
    if ~solfound
        warning("Iteration stopped at after "+num2str(i-1)+...
                               " steps due to numerical problems.")
        return
    end
end
% final analysis:
options.M1M2free = true;
options.rho_iter = 6;
options.rho_start = info(i).rho;
info_ana = iqc_analysis(iqc,G,K,CH,options);
info(i+1).tau = 1;
info(i+1).gam = 0;
if options.ind_ana_each_step
    info(end).ind_ana = info_ana;
end
for j=1:length(CH)
    if info_ana(j).info == "solution found"
        info(i+1).gam = info(i+1).gam + info_ana(j).gam;
        info(i+1).rho(j) = info_ana(j).rho;
    else
        info(i+1).gam = info(i+1).gam + syn_init.gam(j);
        info(i+1).tau = info(i+1).tau_s;
        info(i+1).rho(j) = info(i).rho;
    end
end
info(i+1).time = toc;
end


function G = scale_Delta_by_tau(G, tau)
    p = G.InputGroup.p;
    G.B(:,p) = G.B(:,p)*tau;
    G.D(:,p) = G.D(:,p)*tau;
end

function [t, t_min, t_max] = update_tau(t, t_min, t_max, j, solved)
    if solved
        t_min = t;
    else
        if t_max <= t_min
            t_min = 0;
        else
            t_max = t;
        end
    end
    if j==1
        t = t_min;
    else
        t = (t_max+t_min)/2;
    end
end

function display_info(info)
    i=length(info);
    disp("Iteration "+num2str(i)+":  tau ="+num2str(info(i).tau)+...
         " gamma="+num2str(info(i).gam)+ ...
         " rho="+num2str(info(i).rho))
end