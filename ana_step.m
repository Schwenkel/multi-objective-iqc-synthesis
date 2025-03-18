function [ana_ws,syn_init] = ana_step(iqc,G,K,CH,ana_ws,options)
% set default options if not specified
if ~isfield(options,"rho_min")    % lower bound on rho for which the IQC
    options.rho_min = 0;          % is still valid for Delta_rho
end
if ~isfield(options,"margin_ana") % margin for strict inequalities, i.e.,
    options.margin_ana = 1e-7;    % LMI > 0   <==>   LMI >= margin_pos
end
if ~isfield(options,"M1M2_free")  % specifies whether the restriction (38)
    options.M1M2_free = true;     % for M1, M2, X1 and X2 is used.
end
if nargout == 2
    no_syn=false; % analysis with common P and M with synthesis preparation
    options.M1M2_free = false;   % M1, M2 must satisfy (38) for synthesis
else  % analysis for single channel CH(1) without synthesis preparation
    no_syn=true; 
end
if ~isfield(options,"rho_iter")  % number of line search steps to minimize
    options.rho_iter = 3;        % rho. In iqc_synthesis, several ana_steps
    if no_syn                    % are made, thus, small rho is sufficient.
        options.rho_iter = 15;   % for iqc_analysis we need to be accurate.
    end
end
no_K = false;
if isempty(K)
    no_K = true;
    GK = G;          % K is empty if we only want to analyze G
else
    nx = length(G.A);
    GK = lft(G,K); % get closed loop
end
if isempty(ana_ws)
    ana_ws = struct();
    ana_ws.rho = options.rho_start;
    ana_ws.gam = ones(1,length(CH))*1e5;  % this is only used as a warm 
    ana_ws.mu = ones(1,length(CH))*1e5;   % start
end
rho = ana_ws.rho;
rho_iter = 0;     % no rho optimization
% we distinguish between CH(i).gam that are optimized and those that are a
% constraint. If CH(i).gam is a constraint that is not yet satisfied, then
% it is a priority to optimize this CH(i).gam.
opt_2nd = [];
opt_1st = false;
for j = 1:length(CH)
    if ~isfield(CH(j),"sigma")
        CH(j).sigma=0.5;
    end
    if  CH(j).gam>0 
        if ana_ws.gam(j) > CH(j).gam
            CH(j).gam=-1;   % not yet feasible, hence optimize it
            opt_1st = true;
        else
            ana_ws.gam(j)=CH(j).gam;
            if CH(j).gain ~="p2p"
                ana_ws.mu(j)=CH(j).gam; % mu = gam if not p2p
            end
        end
    else
        opt_2nd = [opt_2nd; j];
    end
    if CH(j).gain == "p2p"  % rho needs to be optimized only for p2p
        rho_iter=options.rho_iter;
        if rho == 1
            rho=0.999;
            warning("rho=1 not possible for p2p, changed to rho=0.999")
        end
    end
end
if opt_1st
    for j = opt_2nd'
        CH(j).gam = 1e4;      % should be larger than the actual value
        ana_ws.gam(j) = 1e4; 
    end
end
eps = options.margin_ana;   % shorthand notation
if isfield(ana_ws,"rhotry")
    rhotry = ana_ws.rhotry;
else
    rhotry = rho-0.01;
end
% Get LMIs and solve the SDP:
[LMIs,margins,ana_var] = get_LMIs(iqc,GK,CH,rho,eps,options.M1M2_free);
assign_warm_start(ana_var, ana_ws);
sdp_opts = sdpsettings('solver','mosek','verbose',0);
optimize(LMIs, sum(ana_var.gam), sdp_opts);
ana_ws = get_values(ana_var, CH(1).sigma);
% Is solution feasible?
error_occured = false;
if check_LMIs(LMIs, margins) < 0
    error_occured=true;
    gam = Inf;
else
    gam = sum(ana_ws.gam);
end
% Optimize rho
for i=1:rho_iter
    [LMIs,margins,ana_var] = get_LMIs(iqc,GK,CH,rhotry,eps,options.M1M2_free);
    assign_warm_start(ana_var, ana_ws);
    optimize(LMIs, sum(ana_var.gam), sdp_opts);
    if check_LMIs(LMIs, margins)<0
        gamtry = Inf;
    else
        error_occured = false;
        ana_ws_try = get_values(ana_var, CH(1).sigma);
        gamtry = sum(ana_ws_try.gam);
        if gamtry < gam  % best gam so far
            ana_ws = ana_ws_try;
        end
    end
    % update rho and find new rhotry
    [rhotry, rho, gam] = step_size_rule(rhotry,gamtry,rho,gam,options.rho_min);
end
if error_occured
    [ana_ws,syn_init] = try_again_with_balanced_K(iqc,G,K,CH,ana_ws,options);
    return
else
    ana_ws.info = "solution found";
end
if no_syn || no_K
    return
end
% Prepare for Synthesis:
[Psih, Xh, Vh, Z] = factorization(ana_ws.Psi1,ana_ws.Psi2,ana_ws.M,ana_ws.X);
iqc_fact = factorized_IQC(Psih, Xh);
% As Psih may have different dimension than Psi, we adjust dimension of K
% as well, as we can only synthesize K of dimension npsi+nx
npsi_1 = size(Psih.A1,1);
npsi_2 = size(Psih.A2,1);
npsih = npsi_1+npsi_2;
K = set_K_dimension(K, npsih+nx);
GK = lft(G,K);
% Find warm start for synthesis
[LMIs, margins, ana_var] = get_LMIs(iqc_fact, GK, CH, rho, eps/2, options.M1M2_free);
optimize(LMIs, sum(ana_var.gam), sdp_opts);
if check_LMIs(LMIs,margins) <0
    [ana_ws,syn_init] = try_again_with_balanced_K(iqc,G,K,CH,ana_ws,options);
    return
end
% struct with synthesis initialization
syn_init = struct();
syn_init.rhotry = rhotry; % remember rhotry for next call of ana_step
syn_init.Psih = Psih;
syn_init.Xh = Xh;
syn_init.H = ana_ws.H;
syn_init.Vh = Vh;
syn_init.Ph = value(ana_var.P);
syn_init.X = ana_ws.X;
syn_init.M = ana_ws.M;
syn_init.Z = Z;
syn_init.K = K;
syn_init.rho = rho;
syn_init.gam = value(ana_var.gam);
syn_init.mu = value(ana_var.mu);
syn_init.X1 = (1-CH(1).sigma)*ana_ws.X;
syn_init.M1 = (1-CH(1).sigma)*ana_ws.M;
syn_init.X2 = CH(1).sigma*ana_ws.X;
syn_init.M2 = CH(1).sigma*ana_ws.M;


    function [ana_ws,syn_init] = try_again_with_balanced_K(iqc,G,K,CH,ana_ws,options)
        syn_init = struct();
        if no_K || K.name == "balanced"
            ana_ws.info = "no solution found";
        else
            K_reduce = reducespec(K,"balanced");
            Kbal = getrom(K_reduce,Order = length(K.A)); % Order may be smaller
            Kbal.name = "balanced";
            ana_ws.P=[];
            if no_syn
                [ana_ws] = ana_step(iqc,G,Kbal,CH,ana_ws,options);
                if ana_ws.info ~= "solution found" && length(Kbal.A)>=1
                    % try again with 1 order reduced
                    Kbal = getrom(K_reduce,Order = length(Kbal.A)-1);
                    Kbal.name = "balanced";
                    ana_ws.P=[];
                    [ana_ws] = ana_step(iqc,G,Kbal,CH,ana_ws,options);
                end
            else
                [ana_ws,syn_init] = ana_step(iqc,G,Kbal,CH,ana_ws,options);
            end
        end
    end

end

% Auxiliary functions:

function [LMIs, margins, ana_var] = get_LMIs(iqc, GK, CH, rho, margin, M1M2_free)
    % loop trafo:
    GK_rho = GK;
    GK_rho.A = rho^-1*GK_rho.A;
    GK_rho.B = rho^-1*GK_rho.B;
    % Building the system Sigma:
    S = build_Sigma(iqc.Psi1,iqc.Psi2,GK_rho);
    % Dimensions and channels
    nP = size(S.A,1);
    p = S.InputGroup.p;
    w = S.InputGroup.w;
    s = S.OutputGroup.s;
    z = S.OutputGroup.z;
    % SDP variables
    if exist("P","var")
        P_ = P;
    else
        P_ = sdpvar(nP);
    end
    gam_ = sdpvar(1, length(CH));
    ana_var = struct();
    ana_var.P = P_;
    ana_var.mu = gam_;  % mu = gam whenever not p2p
    ana_var.H = iqc.H;
    ana_var.Psi1 = iqc.Psi1;
    ana_var.Psi2 = iqc.Psi2;
    ana_var.rho = rho;
    % SDP constraints
    LMIs = iqc.MX_con{1};
    margins = iqc.MX_con{2};
    for j = 1:length(CH)
        Sj = S([s z(CH(j).z)], [p w(CH(j).w)]);
        if CH(j).gam > 0  % feasibility problem with fixed gam
            gam_(j) = CH(j).gam;
        end
        if CH(j).gain == "e2e"
            margins = [margins; margin];
            LMIs = [LMIs, e2e_LMI1(Sj, P_, iqc.M, gam_(j), margins(end))];
        else
            if CH(j).gain == "p2p"
                alpha = rho^2/(1-rho^2);
                ana_var.mu(j) = sdpvar(1); % mu must not be =gam for p2p
                beta = ana_var.mu(j);
                LMIs = [LMIs, gam_(j)>=ana_var.mu(j), ana_var.mu(j)>=0];
                margins = [margins; 0; 0];
            elseif CH(j).gain == "e2p"
                alpha = 1;
                beta = 0;
                ana_var.mu(j) = gam_(j);
            else
                error("Unknown gain specified: "+CH(j).gain)
            end
            if M1M2_free 
                X1j = iqc.X;
                M1j = iqc.M;
                [M2j, X2j, MX2_con, H2] = iqc.get_MX2;
                LMIs = [LMIs, MX2_con{1}];
                margins = [margins; MX2_con{2}];
                ana_var.H = [iqc.H(:); H2(:)];
            else
                X1j = iqc.X*(1-CH(j).sigma);
                X2j = iqc.X*CH(j).sigma;
                M1j = iqc.M*(1-CH(j).sigma);
                M2j = iqc.M*CH(j).sigma;
            end
            margins = [margins; margin; 0.01*margin];
            LMI1 = pe2p_LMI1(Sj, P_, M1j+M2j, ana_var.mu(j), margins(end-1));
            LMI2 = pe2p_LMI2(alpha, beta, Sj, P_, X1j, X2j, M2j, gam_(j), margins(end));
            LMIs = [LMIs, LMI1, LMI2];
        end
    end
    if M1M2_free && exist("M1j","var")
        ana_var.M1 = M1j;
        ana_var.X1 = X1j;
        ana_var.M2 = M2j;
        ana_var.X2 = X2j;
        ana_var.X = X1j+X2j;
    else
        ana_var.M = iqc.M;
        ana_var.X = iqc.X;
    end
    margins = [margins; 10*margin];
    LMIs = [LMIs, PX_pos_LMI(P_, ana_var.X, margins(end))];
    ana_var.gam = gam_;
end

function LMI = e2e_LMI1(S, P, M, gam, margin)  % LMI (20) in 
    ns = length(M);                            % [Schwenkel et. al., 2025]  
    np = length(S.InputGroup.p);
    s = 1:ns;
    zi = ns+1:size(S.C,1);  nzi = length(zi);
    nwi = size(S.B,2)-np;
    nP = size(P,1);
    outer = [eye(nP)          zeros(nP, np+nwi);
             S.A              S.B
             S.C(s,:)         S.D(s,:)
             zeros(nwi, nP)   zeros(nwi,np)  eye(nwi)];
    scaling = 1;%e-2;
    inner = blkdiag(-P, P, M, -gam*scaling*eye(nwi));
    LM = [outer'*inner*outer      [S.C(zi,:)'; S.D(zi,:)'];
          S.C(zi,:)  S.D(zi,:)    -gam/scaling*eye(nzi)];
    LMI = LM <= -margin*eye(length(LM));
end

function LMI = PX_pos_LMI(P, X, margin)       % LMI (15) in
    nP = length(P);                           % [Schwenkel et. al., 2025]  
    LMI = P - blkdiag(X,zeros(nP-length(X))) >= margin*eye(nP);
end

function LMI = pe2p_LMI1(S, P, M, mu, margin)
    ns = length(M);
    np = length(S.InputGroup.p);
    s = 1:ns;
    nwi = size(S.B,2)-np;
    nP = size(P,1);
    outer = [eye(nP)          zeros(nP, np+nwi);
             S.A              S.B
             S.C(s,:)         S.D(s,:)
             zeros(nwi, nP)   zeros(nwi,np)  eye(nwi)];
    inner = blkdiag(-P, P, M, -mu*eye(nwi));
    LM = outer'*inner*outer;
    LMI = LM <= -margin*eye(length(LM));
end

function LMI = pe2p_LMI2(alp, bet, S, P, X1, X2, M2, gam, margin) 
    ns = length(M2);               % LMI (21) in [Schwenkel et. al., 2025]  
    nP = length(P);
    s = 1:ns;
    np = length(S.InputGroup.p);
    zi = ns+1:size(S.C,1);  nzi = length(zi);
    nwi = size(S.B,2)-np;
    outer = [eye(nP)          zeros(nP,np+nwi)
             S.A              S.B
             S.C(s,:)         S.D(s,:)
             zeros(nwi,nP+np)  eye(nwi)];
    X1bar = blkdiag(X1,zeros(nP-length(X1)));
    X2bar = blkdiag(X2,zeros(nP-length(X1)));
    inner = blkdiag(X1bar-P, X2bar, M2, -alp*(gam-bet)*eye(nwi));
    % Schur complement
    LM = [outer'*inner*outer      [S.C(zi,:)'; S.D(zi,:)'];
          S.C(zi,:)  S.D(zi,:)    -gam/alp*eye(nzi)];
    LM = LM + LM'; % numeric noise can cause LM to be non-symmetric
    LMI = LM <= -2*margin*eye(length(LM));
end

function assign_warm_start(ana_var, ana_ws)
    if isfield(ana_ws, "P")
        if length(ana_ws.P) == length(ana_var.P)
            assign(ana_var.P,   ana_ws.P)
        end
        assign(ana_var.mu,  ana_ws.mu)
        assign(ana_var.gam, ana_ws.gam)
        if isfield (ana_var, "M") 
            if ~isfloat(ana_var.X) % If float and no sdpvar "assign" produces an error
                assign(ana_var.X, ana_ws.X)
            end
            if ~isfloat(ana_var.M) % If float and no sdpvar "assign" produces an error
                assign(ana_var.M, ana_ws.M)
            end
            if ~isempty(ana_var.H)
                assign(ana_var.H(:), ana_ws.H(1:length(ana_var.H)));
            end
        end
    end
end

function ana_ws = get_values(ana_var, sigma)
    ana_ws = struct();
    if ~isempty(ana_var.H)
        ana_ws.H = value(ana_var.H);
    else
        ana_ws.H = [];
    end
    ana_ws.P = value(ana_var.P);
    if isfield(ana_var, "M")
        ana_ws.M = value(ana_var.M);
        ana_ws.X = value(ana_var.X);
        ana_ws.M1 = (1-sigma(1))*ana_ws.M;
        ana_ws.X1 = (1-sigma(1))*ana_ws.X;
        ana_ws.M2 = sigma(1)*ana_ws.M;
        ana_ws.X2 = sigma(1)*ana_ws.X;
    else
        ana_ws.M1 = value(ana_var.M1);
        ana_ws.X1 = value(ana_var.X1);
        ana_ws.M2 = value(ana_var.M2);
        ana_ws.X2 = value(ana_var.X2);
        ana_ws.M = ana_ws.M1+ana_ws.M2;
        ana_ws.X = ana_ws.X1+ana_ws.X2;
    end
    ana_ws.mu = value(ana_var.mu);
    ana_ws.gam = value(ana_var.gam);
    ana_ws.Psi1 = ana_var.Psi1;
    ana_ws.Psi2 = ana_var.Psi2;
    ana_ws.rho = ana_var.rho;
end

function feasibility = check_LMIs(LMIs, margins)
    % if feasibility > 0 then feasible
    tol = 1e-7 ;
    feasibility = min(check(LMIs)+margins*4/5+tol);
end

function iqc_fact = factorized_IQC(Psih, Xh)
    nq = size(Psih.B1,2);
    np = size(Psih.B2,2);
    M = blkdiag(eye(nq),-eye(np));
    Psi1 = [ss(Psih.A1,Psih.B1,Psih.C1,Psih.D1,-1); zeros(np,nq)];
    Psi2 = ss(Psih.A2,Psih.B2,[Psih.C3;Psih.C2],[Psih.D3;Psih.D2],-1);
    iqc_fact = IQC(Psi1,Psi2,M,Xh,{[], []},[]);
end

function K = set_K_dimension(Kold, nKnew)
    [nu, ny] = size(Kold.D);
    nKold = length(Kold.A);
    if nKold < nKnew % blow up dimension with unobservable/uncontrollable states
        KoldA = blkdiag(Kold.A, zeros(nKnew-nKold));
        KoldB = [Kold.B; zeros(nKnew-nKold, ny)];
        KoldC = [Kold.C  zeros(nu, nKnew-nKold)];
        K = ss(KoldA, KoldB, KoldC, Kold.D,-1);
    elseif nKold > nKnew
        warning("Controller dimension was reduced")
        K_reduce = reducespec(Kold,"balanced");
        K_rom = getrom(K_reduce,Order = nKnew); % Order may be smaller than nK
        nKrom = size(K_rom.A,1);
        KoldA = blkdiag(K_rom.A, zeros(nKnew-nKrom));
        KoldB = [K_rom.B; zeros(nKnew-nKrom, ny)];
        KoldC = [K_rom.C  zeros(nu, nKnew-nKrom)];
        K = ss(KoldA, KoldB, KoldC, K_rom.D,-1);
    else
        K = Kold;
    end
    K.name = Kold.name;
end

function S = build_Sigma(Psi1,Psi2,GK)
    % channels and dimensions:
    p = GK.InputGroup.p;
    w = GK.InputGroup.w;
    q = GK.OutputGroup.q;
    z = GK.OutputGroup.z;
    nx = length(GK.A);
    np = length(p);
    nw = length(w);
    nz = length(z);
    npsi1 = length(Psi1.A);
    npsi2 = length(Psi2.A);
    ns = size(Psi1.C,1);
    % % Building the system Sigma:
    A = [ Psi1.A              zeros(npsi1,npsi2)   Psi1.B*GK.C(q,:);
          zeros(npsi2,npsi1)  Psi2.A               zeros(npsi2,nx);
          zeros(nx,npsi1)     zeros(nx,npsi2)      GK.A ];
    C = [ Psi1.C             Psi2.C               Psi1.D*GK.C(q,:);
          zeros(nz,npsi1)    zeros(nz,npsi2)      GK.C(z,:) ];
    J = [eye(np) zeros(np,nw)];
    B = [ Psi1.B*GK.D(q,[p w]);
          Psi2.B*J;
          GK.B(:,[p w])];
    D = [ Psi1.D*GK.D(q,[p w])+Psi2.D*J;
          GK.D(z,[p w]); ];
    S = ss(A,B,C,D,-1);
    S.InputGroup.p = p;
    S.InputGroup.w = w;
    S.OutputGroup.s = 1:ns;
    S.OutputGroup.z = ns+1:ns+nz;
end

function [rhotry, rhonew, gam] = step_size_rule(rhotry,gamtry,rho,gam,rhomin)
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
    elseif rhotry<=rhomin
        rhotry = rhomin/2 + rhonew/2;
    end
end
