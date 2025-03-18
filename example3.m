% Example 3:
%% Define model
% parameters:
Ls = 1;             % Shaft lenght
ds = 0.02;          % Shaft diameter
Js = 0;             % Shaft inertia (negligible)
Jm = 0.5;           % Motor inertia
Bm = 0.1;           % Motor viscous friction coefficient
R = 20;             % Resistance of the armature
kt = 10;            % Motor constant
r = 20;             % Gear ratio
k_th = 1280.2;      % Torsional rigidity
Jlinvmax = 1/(10*Jm);      % max inverted load inertia
Jlinvmin = 1/(30*Jm);      % min inverted load inertia
deltamid = (Jlinvmax+Jlinvmin)/2;
deltamax = (Jlinvmax-Jlinvmin)/2;
Bl = 25;            % Load viscous friction coefficient
Vmax =  220; 
Tmax =  78.5398;
% define model
Ts=0.1;
Ac0 = [ 0              1             0                 0; 
       -k_th*deltamid  -Bl*deltamid  k_th/r*deltamid   0;
       0               0             0                 1;
       k_th/(r*Jm)     0             -k_th/(r^2*Jm)    -(Bm + kt^2/R)/Jm ];
Acp = [     0    0       0   0; 
        -k_th  -Bl  k_th/r   0;
            0    0       0   0;
            0    0       0   0 ]; % Ac = Ac0 + delta*Acp;
Bc0 = [0;0;0;kt/(R*Jm)*Vmax];
Bc = [Bc0 zeros(4,2) Bc0];
%% Discretization
% we ignore terms of order delta^2 or higher and of order Ts^4 or higher.
AB0 = expm([Ac0 Bc; zeros(4,8)]*Ts);
Ad0 = AB0(1:4,1:4); Bd0 = AB0(1:4,5:8);
Ad1 = 0*Ac0;
Bd1 = 0*Bd0;
for k = 0:10
    Adj = 0*Ac0;
    for j = 0:k
        Adj = Adj + Ac0^j*Acp*Ac0^(k-j);
    end
    Ad1 = Ad1 + Ts^(k+1)/factorial(k+1)*Adj;
    Bd1 = Bd1 + Ts^(k+2)/factorial(k+2)*Adj*Bc;
end
% Ad = Ad0 + delta*Ad1 + delta^2...;
% Bd0 = Ts*Bc + Ts^2/2*Ac0*Bc + Ts^3/6*Ac0^2*Bc;
% Bd1 = Ts^2/2*Acp*Bc+Ts^3/6*(Acp*Ac0+Ac0*Acp)*Bc;
% Bd = Bd0 + delta*Bd1 + delta^2...;
Cdq1 = Ad1;
Ddq1wup2 = Bd1;
Bdp1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
B = [ Bdp1 Bd0 ];
C = [ Cdq1;    % q1 = J_L*p1 where J_L is unknown
      0          0    0             0;    % q2 = u
      1          0    0             0;    % z1 = x1-w1   tracking error
      k_th/Tmax  0    -k_th/r/Tmax  0;    % z2 = T       torsional torque
      0          0    0             0;    % z3 = u       input voltage
      0          0    0             0;    % y1 = w1      reference
      1          0    0             0 ];  % y = x1+w2    measured output
p1 = 1:4;
p2 = 5;
p = [p1 p2]; % input p
w = 6:7; % input w
u = 8;   % input u
q1 = 1:4;
q2 = 5;
q = [q1 q2]; % output q
z = 6:8; % output z
y = 9:10; % output y
D = zeros(y(end),u(end));
D(z,w) = [-1 0; 0 0; 0 0];
D(y,w) = [1 0; 0 1];
D(q1,[p2 w u]) = Ddq1wup2;
D(q2,u) = 1;
D(z,u) = [0; 0; 1];
%% Performance specifications
CH=struct();
% performance channel 1
CH(1).gam = -1;
CH(1).gain = "e2e";
CH(1).z = 1;
CH(1).w = 1:2;
CH(1).sigma = 0.65;
% performance channel 2
CH(2).gam = -1;
CH(2).gain = "e2p";
CH(2).z = 2;
CH(2).w = 1:length(w);
CH(2).sigma = 0.65;
% performance channel 3
CH(3).gam = -1;
CH(3).gain = "e2p";
CH(3).z = 3;
CH(3).w = 1:length(w);
CH(3).sigma = 0.65;
%% Performance weights:
G = ss(Ad0,B,C,D,Ts);
W1 = makeweight(200,[0.2, 2], 1, Ts)*makeweight(2.4, [3, 1.4], 1, Ts);
W2 = makeweight(1, [1, 4.5], 200, Ts);
W = blkdiag(W1, W2);
np = length(p);
G=G*blkdiag(eye(np),W,eye(length(u)));
G = balreal(G);
G.InputGroup.p = p;
G.InputGroup.w = w;
G.InputGroup.u = u;
G.OutputGroup.q = q;
G.OutputGroup.z = z;
G.OutputGroup.y = y;
%% Uncertainty and IQC
% unmodelled motor dynamics are l2rho bounded:
nu=2;
rhofilter=0.4;
gam_l2 = 0.1;
P0 = [gam_l2 0; 0 -1/gam_l2];
iqc1 = IQC.parametric_interval(nu,rhofilter,-deltamax,deltamax,length(q1));
iqc2 = IQC.dynamic(nu,rhofilter,P0);
iqc = iqc1.combine_with(iqc2);

%% Solve:
opts.rho_start=1;
opts.iterations=11;
opts.ind_ana_each_step = true;
[K_nom, ana_ws] = nom_synthesis(G,CH);
disp("Nominal synthesis: gam="+num2str(ana_ws.gam))
[K,ana_ws,info] = iqc_synthesis(G,CH,iqc,opts);
%% Get real systems
color='k';
%K=K_nom;
zs = tf('z');
pole = 0.4;
DD1 = gam_l2*(1-abs(pole))/(zs+pole);
DD2 = gam_l2*1/(zs);
Delta = {};
Delta{1} = blkdiag(-deltamax*eye(4), DD1);
Delta{2} = blkdiag(-0.5*deltamax*eye(4), -DD1);
Delta{3} = blkdiag(deltamax*eye(4), -DD2);
Delta{4} = blkdiag(-deltamax*eye(4), DD2);
G_no_W = ss(Ad0,B,C,D,Ts);
GKs = {};
GKnoms = {};
for i=1:length(Delta)
    GKs{i} = lft(Delta{i}, lft(G_no_W,K));
    GKnoms{i} = lft(Delta{i}, lft(G_no_W,K_nom));
end
GKs{length(Delta)+1} = lft(G_no_W([z,y],[w,u]),K);
GKnoms{length(Delta)+1} = lft(G_no_W([z,y],[w,u]),K_nom);
%% step response with robust and nominal controller
for i=1:2
    GKis = GKs;
    if i==2
        GKis = GKnoms;
    end
    figure(2+i)
    T = 150;
    rng(1)
    zmax = zeros(3,length(GKis));
    max_meas_noise = 2/180*pi;
    for j=1:length(GKis)
        w2 = randn(1,T)*max_meas_noise;
        w1 = 90/180*pi*ones(1,T);
        GKi = GKis{j};
        xt=zeros(length(GKi.A),T+1);
        zt=zeros(size(GKi.C,1),T);
        for t=1:T
            xt(:,t+1) = GKi.A*xt(:,t)+GKi.B*[w1(t); w2(t)];
            zt(:,t) = GKi.C*xt(:,t)+GKi.D*[w1(t); w2(t)];
            zt(1,t) = (zt(1,t) + w1(t))/pi*180; % output instead of error
        end
        zt(3,:) = zt(3,:);
        subplot 131
        stairs((1:T)*Ts, zt(1,:),color);hold on
        subplot 132
        stairs((1:T)*Ts, zt(2,:),color);hold on
        subplot 133
        stairs((1:T)*Ts, zt(3,:),color);hold on
        zmax(1,j) = max(abs(zt(1,:)));
        zmax(2,j) = max(abs(zt(2,:)));
        zmax(3,j) = max(abs(zt(3,:)));
    end
end
%% bode plots
figure(1)
 bode(1/W1,'r')
hold on;
for i = 1:length(GKs)
    bode(GKs{i}(1,1),color)
end
figure(2)
 bode(1/W2,'r')
hold on;
for i = 1:length(GKs)
    bode(GKs{i}(1,2),color)
end
%% gamma over iterations
figure(5)
gam_a = info.gam;
for i=2:length(info)
    if ~isempty(info(i).gam)
        plot(i,info(i).gam,'kx')
    end
    hold on
end
for i=2:length(info)
    ind_gam = 0;
    for j=1:length(CH)
        ind_gam=ind_gam+info(i).ind_ana(j).gam;
    end
    plot(i,ind_gam,'rx')
end
%% average computation time
time = 0;
for i=1:length(info)
    time = time + info(i).time;
end
disp("Average computation time per iteration: "+(time/length(info)));
% note that the individual analysis is not added to the computation time
% of each iteration as it is not needed for the iteration and just computed
% for comparison.
%% Comparison of reduced order controllers
figure(6)
K_reduce = reducespec(K,"balanced");
gammas = zeros(length(CH), length(K.A));
for i = 0:length(K.A)-1
    K_rom = getrom(K_reduce,Order = i);
    infoi = iqc_analysis(iqc,G,K_rom,CH,opts);
    for j = 1:length(CH)
        if infoi(j).info == "solution found"
            gammas(j,i+1) = infoi(j).gam;
        else
            gammas(j,i+1) = inf;
        end
    end
end
plot(sum(gammas))
ylim([3.5,4.5])