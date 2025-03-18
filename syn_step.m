function [K, ana_ws] = syn_step(G,CH,syn_init,options)
if ~isfield(options,"margin_syn") % margin for strict inequalities, i.e.,
    options.margin_syn = 5e-8;    % LMI > 0   <==>   LMI >= margin_pos
end
eps = options.margin_syn;
options = sdpsettings('solver','mosek','verbose',0);
ana_ws = struct();

% dimensions
nu=length(G.InputGroup.u);
ny=length(G.OutputGroup.y);
nx = size(G.A,1);
n1 = size(syn_init.Psih.A1,1);
n2 = size(syn_init.Psih.A2,1);
npsih = n1+n2;
nK = nx+npsih;

rho = syn_init.rho;
Grho = G;
Grho.A = G.A/rho;
Grho.B = G.B/rho;
Theta = construct_generalized_plant(Grho,syn_init.Psih);

% LMI decision variables for synthesis
syn_var.P22 = sdpvar(nK);
syn_var.P11 = sdpvar(nK);
syn_var.KLMN = sdpvar(nK+nu, nK+ny, 'full');
syn_var.gam = sdpvar(1,length(CH));
syn_var.mu = syn_var.gam; % mu = gam for all but p2p
opt_2nd = [];
opt_1st = false;
for i = 1:length(CH)
    if CH(i).gam > 0 
        if syn_init.gam(i) > CH(i).gam
            CH(i).gam=-1; % not yet feasible, hence optimize it
            opt_1st = true;
        else
            syn_var.gam(i) = CH(i).gam;% feasibility problem with fixed gam
            syn_var.mu(i) = CH(i).gam;
            syn_init.gam(i) = CH(i).gam;
        end
    else
        opt_2nd = [opt_2nd; i];
    end
    if CH(i).gain == "p2p"
        syn_var.mu(i)=sdpvar(1); % mu not equal gam for p2p
        if syn_init.rho == 1
            syn_init.rho=0.999;
            warning("rho=1 not possible for p2p, changed to rho=0.999")
        end
    end
end
if opt_1st
    for i = opt_2nd'
        CH(i).gam = 1e3;
        syn_var.gam(i) = CH(i).gam;
        if CH(i).gain ~= "p2p"
            syn_var.mu(i) = CH(i).gam;
        end
        syn_init.gam(i) = CH(i).gam;
    end
end

% warm starts for the optimization
if CH(1).gain == "init"  % in the initialization there is no warm start K
    syn_ws=struct();
    syn_ws.rho = syn_init.rho;
else
    syn_ws = assign_warm_start(syn_var,Theta,syn_init);
end

[LMIs, margins] = get_LMIs(Theta, CH, syn_var, syn_ws, eps);

if check_LMIs(LMIs,margins)<0%-1e-5
    disp("Warning: Warmstart for synthesis not feasible by " ...
                   +num2str(-check_LMIs(LMIs,margins))+".")
end
sol = optimize(LMIs,sum(syn_var.gam), options);
if check_LMIs(LMIs,margins)<0
    ana_ws.info = "no solution found";
    K=[];
    return
end
syn_ws = get_values(syn_var);

% Numerical improvevements from Remark 4.6 in [Scherer, Weiland, 2005]
%
% fix almost optimal gam
gam = 1.001*sum(syn_ws.gam);
% step 1: small X, Y, K, L, M, N
a1 = sdpvar(1);
LMIs = [LMIs, sum(syn_var.gam)<= gam, bound_by_a1(a1,syn_var)];
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
LMIa2 = [syn_var.P11 a2*eye(nK); a2*eye(nK) syn_var.P22] >= margins(end)*eye(2*nK);
LMIs = [LMIs, a1<=a1max, LMIa2];
optimize(LMIs, -a2, options);
if check_LMIs(LMIs, margins) < 0
    ana_ws.info = "no solution found";
    K = [];
    return
else
    ana_ws.info = "solution found";
    % Obtain controller matrices K = (A,B,C,D)
    syn_ws = get_values(syn_var);
end

[Krho, Ph] = obtain_controller(Theta, syn_ws);
K = Krho;
K.A = K.A*rho;
K.B = K.B*rho;
ana_ws.gam = syn_ws.gam;
ana_ws.mu = syn_ws.mu;
ana_ws.M = syn_init.M;
ana_ws.X = syn_init.X;
ana_ws.H = syn_init.H;
ana_ws.rho = syn_init.rho;
ana_ws.rhotry = syn_init.rhotry;
ana_ws.P = Ph2P(Ph,syn_init.Vh,syn_init.Z);
end

% Auxiliary functions


function [LMIs, margins] = get_LMIs(Theta,CH,syn_var,syn_ws, margin)
    nK = length(Theta.A);
    CL = get_syn_CL(Theta,syn_var);
    CL_old = get_syn_CL(Theta,syn_ws);
    P_ = [ syn_var.P11    eye(nK)
           eye(nK)        syn_var.P22    ];
    P11o = syn_ws.P11; % old P11
    Xh = syn_ws.Xh;
    % SDP constraints
    margins = margin;
    LMIs = PX_pos_LMI(P_,P11o,Xh,margins(end));
    for j = 1:length(CH)
        CLj = get_syn_CLj(CL,CH(j));
        gamj = syn_var.gam(j);
        if CH(j).gain == "e2e"
            margins(end+1,1) = margin;
            LMI2 = e2e_LMI1(P_, CLj, gamj, margins(end));
            LMIs = [LMIs, LMI2];
        else
            if CH(j).gain == "p2p"
                a = syn_ws.rho^2/(1-syn_ws.rho^2);
                b = syn_var.mu(j);
                LMIs = [LMIs, syn_var.gam(j) >= syn_var.mu(j), syn_var.mu(j)>=0];
                margins = [margins;0;0];
            elseif CH(j).gain == "e2p"
                a = 1;
                b = 0;
            else
                error("Unknown gain specified: "+CH(j).gain)
            end
            CLjo = get_syn_CLj(CL_old,CH(j));
            sj  = CH(j).sigma;
            margins = [margins; eps; 0.01*eps];
            LMI2 = pe2p_LMI1(P_,CLj,syn_var.mu(j),margins(end-1));
            LMI3 = pe2p_LMI2(a,b,P_,P11o,Xh,CLj,CLjo,sj,gamj,margins(end));
            LMIs = [LMIs, LMI2, LMI3];
        end
    end
end

function LMI = pe2p_LMI1(P, CL, mu, margin)
    % according to (53) in [Schwenkel et. al., 2025]  
    np = length(CL.InputGroup.p);
    nw = size(CL.B,2)-np;
    s = CL.OutputGroup.s;   ns =  length(s);
    LMI11 = blkdiag(-P, -eye(np), -mu*eye(nw));
    LM21 = [CL.A  CL.B; CL.C(s,:) CL.D(s,:)];
    LMI22 = blkdiag(-P, -eye(ns));
    LMI = [ LMI11  LM21'; 
            LM21  LMI22 ] <= - margin*eye(length(LMI11)+length(LMI22));
end

function LMI = e2e_LMI1(P, CL, gam, margin)
    % according to (52) in [Schwenkel et. al., 2025]  
    np = length(CL.InputGroup.p);
    nw = size(CL.B,2)-np;
    ns = length(CL.OutputGroup.s);
    nz = size(CL.C,1)-ns;
    LMI11 = blkdiag(-P, -eye(np), -gam*eye(nw));
    LM21 = [CL.A CL.B; CL.C, CL.D];
    LMI22 = blkdiag(-P, -eye(ns), -gam*eye(nz));
    LMI = [ LMI11  LM21'; 
            LM21   LMI22 ] <= - margin*eye(length(LMI11)+length(LMI22));
end

function LMI = PX_pos_LMI(P, P11old, Xh, margin)
    % according to (51) in [Schwenkel et. al., 2025]  
    nP = length(P);
    nK = nP/2;
    P11 = P(1:nK, 1:nK);
    [L, X2] = decomposeZ(Xh);
    [nL, npsih] = size(L);
    E11 = [P11 eye(nK)];
    E11 = E11(1:npsih, :);
    E11_old = [P11old eye(nK)];
    E11_old = E11_old(1:npsih, :);
    EXE11 = E11_old'*X2*E11_old-2*E11_old'*X2*E11;
    EXE11 = (EXE11+EXE11')/2;
    LMI = [P-EXE11    E11'*L'; 
           L*E11      eye(nL) ] >= margin*eye(nP+nL);
end

function LMI = pe2p_LMI2(a,b,P,P11o,Xh,CL,CLo,sigma,gam,margin)
    % according to (54) in [Schwenkel et. al., 2025]  
    ns = length(CL.OutputGroup.s);
    zj = ns+1:size(CL.C,1); nzj = length(zj);
    np = length(CL.InputGroup.p);
    nwj = size(CL.B,2)-np;
    nK = length(P)/2;
    P11 = P(1:nK, 1:nK);
    [L, X2] = decomposeZ(Xh);
    [nL, npsih] = size(L);
    E11 = [P11 eye(nK)];
    E11 = E11(1:npsih, :);
    E11_old = [P11o eye(nK)];
    E11_old = E11_old(1:npsih, :);
    E1 = [E11 zeros(npsih, np+nwj)];
    E2 = [CL.A(1:npsih, :) CL.B(1:npsih, :)];
    E1_old = [E11_old zeros(npsih, np+nwj)];
    E2_old = [CLo.A(1:npsih, :) CLo.B(1:npsih, :)];
    EXE1 = E1_old'*X2*E1_old-2*E1_old'*X2*E1;
    EXE1 = (EXE1+EXE1')/2;
    EXE2 = E2_old'*X2*E2_old-2*E2_old'*X2*E2;
    EXE2 = (EXE2+EXE2')/2;
    % LMI for sigma=0 and sigma=1 differs as discussed in Remark 9 in 
    % [Schwenkel et. al., 2025].
    sigma_tol=1e-5;
    if (1-sigma) < sigma_tol % i.e. sigma==0
        LMI11 = blkdiag(-P, -sigma*eye(np), -a*(gam-b)*eye(nwj));
        LMI11 = LMI11 + EXE2;
        LMI21 = [L*E2; CL.C  CL.D];
        LMI22 = blkdiag(-1/sigma*eye(nL),-1/sigma*eye(ns),-gam/a*eye(nzj));
    elseif sigma < sigma_tol % i.e. sigma==1
        T = [ eye(2*nK) zeros(2*nK,nwj); 
              zeros(np, 2*nK+nwj); 
              zeros(nwj,2*nK) eye(nwj) ];
        LMI11 = blkdiag(-P, -a*(gam-b)*eye(nwj));
        LMI11 = LMI11 + T'*EXE1*T;
        LMI21 = [L*E1; CL.C(zj,:)  CL.D(zj,:)]*T;
        LMI22 = blkdiag(-1/(1-sigma)*eye(nL), -gam/a*eye(nzj));
    else
        LMI11 = blkdiag(-P, -sigma*eye(np), -a*(gam-b)*eye(nwj));
        LMI11 = LMI11 + (1-sigma)*EXE1 + sigma*EXE2;
        LMI21 = [L*E1; L*E2; CL.C  CL.D];
        LMI22 = blkdiag(-1/(1-sigma)*eye(nL), -1/sigma*eye(nL), ...
                        -1/sigma*eye(ns), -gam/a*eye(nzj));
    end
    nLMI = length(LMI11)+length(LMI22);
    LMI = [LMI11 LMI21'; LMI21 LMI22] <= -margin*eye(nLMI);
end


function syn_CL = get_syn_CL(Ggen,syn_param)
    % according to (49) in [Schwenkel et. al., 2025]  
    % syn_param can be syn_var or syn_ws
    nK = length(Ggen.A);
    p = Ggen.InputGroup.p;
    w = Ggen.InputGroup.w;
    u = Ggen.InputGroup.u;
    s = Ggen.OutputGroup.s;
    z = Ggen.OutputGroup.z;
    y = Ggen.OutputGroup.y;
    P11 = syn_param.P11;
    P22 = syn_param.P22;
    A = Ggen.A;
    Bu = Ggen.B(:,u);
    Bpw = Ggen.B(:,[p w]);
    Cy = Ggen.C(y,:);
    Csz = Ggen.C([s z],:);
    Dszu = Ggen.D([s z],u);
    Dszpw = Ggen.D([s z],[p w]);
    Dypw = Ggen.D(y,[p w]);
    KLMN = syn_param.KLMN;
    K = KLMN(1:nK,1:nK);
    M = KLMN(nK+1:end, 1:nK);
    L = KLMN(1:nK,nK+1:end);
    N = KLMN(nK+1:end,nK+1:end);
    Anu = [A*P11+Bu*M    A+Bu*N*Cy 
           K             P22*A+L*Cy   ];
    Bnu = [Bpw+Bu*N*Dypw
           P22*Bpw+L*Dypw   ];
    Cnu = [Csz*P11+Dszu*M   Csz+Dszu*N*Cy];
    Dnu = Dszpw+Dszu*N*Dypw;
    syn_CL = struct(); % ss model is not possible if KLMN is sdpvar
    syn_CL.A = Anu; syn_CL.B = Bnu; syn_CL.C = Cnu; syn_CL.D = Dnu;
    syn_CL.InputGroup = Ggen.InputGroup;
    syn_CL.OutputGroup = Ggen.OutputGroup;
end

function LMIs = bound_by_a1(a1,syn_var)
    nK = length(syn_var.P22);
    [nKu, nKy] = size(syn_var.KLMN);
    LMIs = syn_var.P22 <= a1*eye(nK);
    LMIs = [LMIs, syn_var.P11 <= a1*eye(nK)];
    LMIs = [LMIs, [a1*eye(nKu)    syn_var.KLMN; 
                   syn_var.KLMN'    a1*eye(nKy)  ] >=0 ];
end
    
function Ggen = construct_generalized_plant(G,Psih)  
    % according to (45) in [Schwenkel et. al., 2025]  
    n1 = length(Psih.A1);
    n2 = length(Psih.A2);
    nx = length(G.A);
    p = G.InputGroup.p;
    w = G.InputGroup.w;    nw = length(w);
    u = G.InputGroup.u;    nu = length(u);
    q = G.OutputGroup.q;
    z = G.OutputGroup.z;   nz = length(z);
    y = G.OutputGroup.y;   ny = length(y);
    C2i = -Psih.D2\Psih.C2;
    A2i = Psih.A2+Psih.B2*C2i;
    B2i = Psih.B2/Psih.D2;
    A = [ Psih.A1        Psih.B1*G.D(q,p)*C2i   Psih.B1*G.C(q,:)
          zeros(n2,n1)   A2i                    zeros(n2,nx)
          zeros(nx,n1)   G.B(:,p)*C2i           G.A          ];
    B = [ Psih.B1*(G.D(q,p)/Psih.D2)  Psih.B1*G.D(q,w)  Psih.B1*G.D(q,u)
          B2i                         zeros(n2,nw)      zeros(n2,nu)
          G.B(:,p)/Psih.D2            G.B(:,w)          G.B(:,u)   ];
    Dsp = (Psih.D1*G.D(q,p)+Psih.D3);
    C = [ Psih.C1       Psih.C3+Dsp*C2i  Psih.D1*G.C(q,:)
          zeros(nz,n1)  G.D(z,p)*C2i     G.C(z,:)
          zeros(ny,n1)  G.D(y,p)*C2i     G.C(y,:) ];
    D = [ Dsp/Psih.D2       Psih.D1*G.D(q,w)  Psih.D1*G.D(q,u)
          G.D(z,p)/Psih.D2  G.D(z,w)          G.D(z,u)  
          G.D(y,p)/Psih.D2  G.D(y,w)          G.D(y,u) ];
    Ggen = ss(A,B,C,D,-1);
    Ggen.InputGroup = G.InputGroup;
    Ggen.OutputGroup.s = G.OutputGroup.q;  
    Ggen.OutputGroup.z = G.OutputGroup.z;  
    Ggen.OutputGroup.y = G.OutputGroup.y;  
end

    function syn_ws = get_controller_syn_ws(Phat, Theta, K)
        u = Theta.InputGroup.u;
        y = Theta.OutputGroup.y;
        nK = size(K.A,1);
        [nu, ny] = size(K.D);
        syn_ws.P22 = Phat(1:nK,1:nK);
        syn_ws.Uold = Phat(1:nK,nK+1:end);
        Phatinv = inv(Phat);
        syn_ws.P11 = Phatinv(1:nK,1:nK);
        syn_ws.Vold = Phatinv(1:nK,nK+1:end);
        UXBI = [syn_ws.Uold  syn_ws.P22*Theta.B(:,u); zeros(nu, nK) eye(nu)];
        VCYI = [syn_ws.Vold' zeros(nK, ny); Theta.C(y,:)*syn_ws.P11  eye(ny)];
        KLMN = UXBI*[K.A K.B; K.C K.D]*VCYI;
        KLMN = KLMN + blkdiag(syn_ws.P22*(Theta.A*syn_ws.P11), zeros(nu, ny));
        syn_ws.KLMN = KLMN;
    end

    function [K, Ph] = obtain_controller(Theta, syn_ws)
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
        Y = [P11 V; eye(nK) zeros(nK)]';
        Ph = Y'\[eye(nK) zeros(nK); P22 U];
    end

    function syn_ws = assign_warm_start(syn_var,Theta,syn_init)
        Krho = syn_init.K;
        Krho.A = Krho.A/syn_init.rho;
        Krho.B = Krho.B/syn_init.rho;
        syn_ws = get_controller_syn_ws(syn_init.Ph, Theta, Krho);
        syn_ws.mu = syn_init.mu;
        syn_ws.gam = syn_init.gam;
        syn_ws.Xh = syn_init.Xh;
        syn_ws.rho = syn_init.rho;
        assign(syn_var.KLMN, syn_ws.KLMN)
        assign(syn_var.P11, syn_ws.P11)
        assign(syn_var.P22, syn_ws.P22)
        assign(syn_var.mu,  syn_ws.mu)
        assign(syn_var.gam, syn_ws.gam)
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
        tol = 5e-8;
        feasibility = min(check(LMIs)+margins*1/2+tol);
    end

function CLj = get_syn_CLj(syn_CL,CHj)
    p = syn_CL.InputGroup.p;
    w = syn_CL.InputGroup.w;
    s = syn_CL.OutputGroup.s;
    z = syn_CL.OutputGroup.z;
    CLj.A = syn_CL.A;
    CLj.B = syn_CL.B(:,[p w(CHj.w)]);
    CLj.C = syn_CL.C([s z(CHj.z)], :);
    CLj.D = syn_CL.D([s z(CHj.z)], [p w(CHj.w)]);
    CLj.InputGroup.p = syn_CL.InputGroup.p;
    CLj.OutputGroup.s = syn_CL.OutputGroup.s;

end


function P = Ph2P(Ph, Vh, Z)
    nPh = length(Ph);
    [npsi, npsih] = size(Vh);
    Vperp = null(Vh);
    T = blkdiag([Vperp'; Vh], eye(nPh-npsih));
    Pbar = T'\((Ph - blkdiag(Z, zeros(nPh-npsih)))/T);
    ndiff = npsih-npsi;
    P = Pbar(ndiff+1:end,ndiff+1:end)-Pbar(ndiff+1:end,1:ndiff)*(Pbar(1:ndiff,1:ndiff)\Pbar(1:ndiff,ndiff+1:end));
    P = (P+P')/2;  % ensure symmetry
end

function [L, X2] = decomposeZ(X)
    X = (X+X')/2; % ensure symmetry
    [V, D] = eig(X);
    tol = 1e-10;
    L = sqrt(D(diag(D)>tol,:))*V'; % X1=L'*L
    X2 = V*(D.*(-(D<-tol)))*V';
    X2 = (X2+X2')/2; % ensure symmetry
    % X=X1-X2
end