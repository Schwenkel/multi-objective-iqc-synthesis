function [K, info] = nom_synthesis(G,CH)
% solves the nominal mixed synthesis problem
eps = 5e-7; % margin for strict inequalities
options = sdpsettings('solver','mosek','verbose',0);
info = struct();
w = G.InputGroup.w;
u = G.InputGroup.u;
z = G.OutputGroup.z;
y = G.OutputGroup.y;
G = G([z y], [w u]);
% dimensions
nu=length(G.InputGroup.u);
ny=length(G.OutputGroup.y);
nx = size(G.A,1);

rho = 1;

% LMI decision variables for synthesis
syn_var.P22 = sdpvar(nx);
syn_var.P11 = sdpvar(nx);
syn_var.KLMN = sdpvar(nx+nu, nx+ny, 'full');
syn_var.gam = sdpvar(1,length(CH));
syn_var.mu = syn_var.gam; % mu = gam for all but p2p
gam_ = 0;
rho_iter=1;
for i = 1:length(CH)
    if CH(i).gam > 0       
        syn_var.gam(i) = CH(i).gam;% feasibility problem with fixed gam
        syn_var.mu(i) = CH(i).gam;
    end
    if CH(i).gain == "p2p"
        syn_var.mu(i)=sdpvar(1); % mu not equal gam for p2p
        rho=0.99;
        rho_iter=15;
    end
    if isfield(CH(i),"weight")
        gam_ = gam_ + CH(i).weight*syn_var.gam(i);
    else
        gam_ = gam_ + syn_var.gam(i);
    end
end
rhotry = rho;
gam = inf;
for i=1:rho_iter
    Grho = G;
    Grho.A = G.A/rhotry;
    Grho.B = G.B/rhotry;
    [LMIs, margins] = get_LMIs(Grho, CH, syn_var, rhotry, eps);
    optimize(LMIs,gam_, options);
    if check_LMIs(LMIs,margins)<0
        gamtry = Inf;        
    else
        syn_ws_try = get_values(syn_var);
        gamtry = value(gam_);
        if gamtry < gam  % best gama so far
            syn_ws = syn_ws_try;
        end
    end
    if i==1
        gam = gamtry;
        rhotry = rhotry - 0.01;
    else
        [rhotry, rho, gam] = step_size_rule(rhotry,gamtry,rho,gam);
    end
end
if gam == inf
    error("No solution found")
end
% Numerical improvevements from Remark 4.6 in [Scherer, Weiland, 2005]
% fix almost optimal gam
gam = 1.01*value(gam_);
if CH(1).gain == "init"
    gam = 10*gam;
end
% step 1: small X, Y, K, L, M, N
a1 = sdpvar(1);
LMIs = [LMIs, gam_<= gam, bound_by_a1(a1,syn_var)];
margins = [margins; 0; 0;0;0];
optimize(LMIs, a1, options);
if check_LMIs(LMIs,margins)<0
    optimize(bound_by_a1(a1,syn_ws),a1,options);
    a1max = 1.1*value(a1);
else
    a1max = 1.1*value(a1);
    syn_ws = get_values(syn_var);
end
% step 2: push I-XY away from being singular
a2 = sdpvar(1);
margins = [margins; 0; eps];
LMIa2 = [syn_var.P11 a2*eye(nx); a2*eye(nx) syn_var.P22] >= margins(end)*eye(2*nx);
LMIs = [LMIs, a1<=a1max, LMIa2];
optimize(LMIs, -a2, options);
if check_LMIs(LMIs, margins) < 0
    if CH(1).gain == "init"
        info.info = "no solution found";
        return
    else
        disp("WARNING: numerical problems occured when recovering the controller parameters.")
    end
else
    info.info = "solution found";
    % Obtain controller matrices K = (A,B,C,D)
    syn_ws = get_values(syn_var);
end

Krho = obtain_controller(Grho, syn_ws);
K = Krho;
K.A = K.A*rho;
K.B = K.B*rho;
info.gam = syn_ws.gam;
info.mu = syn_ws.mu;
info.sigma = 0.5;
info.rho = rho;

end

% Auxiliary functions


function [LMIs, margins] = get_LMIs(Grho,CH,syn_var,rho, margin)
    nx = length(Grho.A);
    CL     = get_syn_CL(Grho,syn_var);
    P_  = [ syn_var.P11    eye(nx)
            eye(nx)        syn_var.P22    ];
    % SDP constraints
    margins=[];LMIs=[];
    for j = 1:length(CH)
        CLj = get_syn_CLj(CL,CH(j));
        gamj = syn_var.gam(j);
        margins(end+1,1) = margin;
        LMI1 = PX_pos_LMI(P_,margins(end));
        if CH(j).gain == "e2e"
            margins(end+1,1) = margin;
            LMI2 = e2e_LMI1(P_, CLj, gamj, margins(end));
            LMIs = [LMIs, LMI1, LMI2];
        else
            if CH(j).gain == "p2p"
                a = rho^2/(1-rho^2);
                b = syn_var.mu(j);
            elseif CH(j).gain == "e2p"
                a = 1;
                b = 0;
            else
                error("Unknown gain specified: "+CH(j).gain)
            end
            margins = [margins; eps; 0.01*eps];
            LMI2 = pe2p_LMI1(P_,CLj,syn_var.mu(j),margins(end-1));
            LMI3 = pe2p_LMI2(a,b,P_,CLj,gamj,margins(end));
            LMIs = [LMIs, LMI1, LMI2, LMI3];
        end
    end
end

function LMI = pe2p_LMI1(P, CLj, muj, margin)
    nw = size(CLj.B,2);
    LMI11 = blkdiag(-P, -muj*eye(nw));
    LM21 = [CLj.A  CLj.B];
    LMI22 = -P;
    LMI = [ LMI11  LM21'; 
            LM21  LMI22 ] <= - margin*eye(length(LMI11)+length(LMI22));
end

function LMI = e2e_LMI1(P, CLj, gamj, margin)
    nw = size(CLj.B,2);
    nz = size(CLj.C,1);
    scaling=1;%e-2;
    LMI11 = blkdiag(-P, -gamj*scaling*eye(nw));
    LM21 = [CLj.A CLj.B; CLj.C, CLj.D];
    LMI22 = blkdiag(-P, -gamj/scaling*eye(nz));
    LMI = [ LMI11  LM21'; 
            LM21   LMI22 ] <= - margin*eye(length(LMI11)+length(LMI22));
end

function LMI = PX_pos_LMI(P, margin)
    LMI = P >= margin*eye(length(P));
end

function LMI = pe2p_LMI2(a,b,P,CLj,gam,margin)
    nzj = size(CLj.C,1);
    nwj = size(CLj.B,2);
    LMI11 = blkdiag(-P, -a*(gam-b)*eye(nwj));
    LMI21 = [CLj.C  CLj.D];
    LMI22 = -gam/a*eye(nzj);
    nLMI = length(LMI11)+length(LMI22);
    LMI = [LMI11 LMI21'; LMI21 LMI22] <= -margin*eye(nLMI);
end


function syn_CL = get_syn_CL(G,syn_param)
    % syn_param can be syn_var or syn_ws
    nx = length(G.A);
    w = G.InputGroup.w;
    u = G.InputGroup.u;
    z = G.OutputGroup.z;
    y = G.OutputGroup.y;
    P11 = syn_param.P11;
    P22 = syn_param.P22;
    A = G.A;
    Bu = G.B(:,u);
    Bw = G.B(:,w);
    Cy = G.C(y,:);
    Cz = G.C(z,:);
    Dzu = G.D(z,u);
    Dzw = G.D(z,w);
    Dyw = G.D(y,w);
    KLMN = syn_param.KLMN;
    K = KLMN(1:nx,1:nx);
    M = KLMN(nx+1:end, 1:nx);
    L = KLMN(1:nx,nx+1:end);
    N = KLMN(nx+1:end,nx+1:end);
    Anu = [A*P11+Bu*M    A+Bu*N*Cy 
           K             P22*A+L*Cy   ];
    Bnu = [Bw+Bu*N*Dyw
           P22*Bw+L*Dyw   ];
    Cnu = [Cz*P11+Dzu*M   Cz+Dzu*N*Cy];
    Dnu = Dzw+Dzu*N*Dyw;
    syn_CL = struct(); % ss model is not possible if KLMN is sdpvar
    syn_CL.A = Anu; syn_CL.B = Bnu; syn_CL.C = Cnu; syn_CL.D = Dnu;
    syn_CL.InputGroup = G.InputGroup;
    syn_CL.OutputGroup = G.OutputGroup;
end

function LMIs = bound_by_a1(a1,syn_var)
    nK = length(syn_var.P22);
    [nKu, nKy] = size(syn_var.KLMN);
    LMIs = syn_var.P22 <= a1*eye(nK);
    LMIs = [LMIs, syn_var.P11 <= a1*eye(nK)];
    LMIs = [LMIs, [a1*eye(nKu)    syn_var.KLMN; 
                   syn_var.KLMN'    a1*eye(nKy)  ] >=0 ];
end

    function K = obtain_controller(Theta, syn_ws)
        u = Theta.InputGroup.u;   nu = length(u);
        y = Theta.OutputGroup.y;  ny = length(y);
        nK = length(Theta.A);
        % Decomposition I-XY = UV'   (4.2.14)
        % any decpomposition works, we use one where U and V have same min sv
        P11 = syn_ws.P11;
        P22 = syn_ws.P22;
        [Ud, Sd, Vd] = svd(eye(nK) - P22*P11);
        Sd_sqrt_diag = sqrt(diag(Sd));
        if min(Sd_sqrt_diag)<1e-4
            warning("Transformation of synthesis to controller " + ...
                    "variables almost singular.")
        end
        U = Ud*diag(Sd_sqrt_diag);
        V = Vd*diag(Sd_sqrt_diag);
        KABCD = syn_ws.KLMN - blkdiag(P22*Theta.A*P11, zeros(nu,ny));
        KABCD = [ U               P22*Theta.B(:,u) 
                  zeros(nu,nK)     eye(nu) ]\KABCD;
        KABCD = KABCD/[ V'             zeros(nK, ny) 
                        Theta.C(y,:)*P11   eye(ny)       ];
        K = ss(KABCD(1:nK, 1:nK),     KABCD(1:nK, nK+1:end), ...
               KABCD(nK+1:end, 1:nK), KABCD(nK+1:end, nK+1:end), -1);
    end

    function syn_ws = get_values(syn_var)
        syn_ws.KLMN = value(syn_var.KLMN);
        syn_ws.P11 = value(syn_var.P11);
        syn_ws.P22 = value(syn_var.P22);
        syn_ws.mu = value(syn_var.mu);
        syn_ws.gam = value(syn_var.gam);
    end

    function feasibility = check_LMIs(LMIs, margins)
        % if feasibility > 0 then feasible
        tol = 1e-7;
        feasibility = min(check(LMIs)+margins*4/5+tol);
    end

function CLj = get_syn_CLj(syn_CL,CHj)
    w = syn_CL.InputGroup.w;
    z = syn_CL.OutputGroup.z;
    CLj.A = syn_CL.A;
    CLj.B = syn_CL.B(:,w(CHj.w));
    CLj.C = syn_CL.C(z(CHj.z), :);
    CLj.D = syn_CL.D(z(CHj.z), w(CHj.w));

end


function [rhotry, rhonew, gam] = step_size_rule(rhotry,gamtry,rho,gam)
    if gamtry<gam  % direction was good, go 2x in that direction
        rhonew = rhotry;
        gam = gamtry;
        rhotry = rhotry+2*(rhotry-rho);
    else % direction was bad, go 0.33x in opposite direction
        rhonew = rho;
        rhotry = rho+(rho-rhotry)/3;
    end
    if abs(rhotry-rhonew) < 1e-4 % make sure steps do not get too small
        rhotry = rhonew+sign(rhotry-rhonew)*5e-2;
    end
    if rhotry>=1 % ensure rho in (rhomin,1)
        rhotry = 1/2+rhonew/2;
    elseif rhotry<=0
        rhotry = rhonew/2;
    end
end
