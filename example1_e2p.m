% Example 1 in [Schwenkel et. al., 2025]
% Revisits Example Section IV. B from 
%
%   Abou Jaoude, and Farhood, Guaranteed output bounds using performance
%     integral quadratic constraints, Proc. IEEE American Control Conf.,
%     2020
% 
% Therein, a value of gam = 2.683 is achieved as upper bound to the
% e2p-gain.
A = -0.5;
B = [0.5 0.4];
C = [2.5; 2];
D = [0 0.6; 0 0.9];
G = ss(A,B,C,D,-1);
G.InputGroup.p = 1;
G.InputGroup.w = 2;
G.OutputGroup.q = 1;
G.OutputGroup.z = 2;
Psi1 = ss(-0.3, 1.3, [0; -0.1; 0; 0], [0.2; 0; -0.5; 0]);
Psi2 = ss([], [], [], [0; -0.1; 0.3; 1.7]);
lambda1 = sdpvar(2,1);
M = diag([lambda1(1), -lambda1(1), lambda1(2), -lambda1(2)]);
X=0;
MX_con = {lambda1>=0, 0};
iqc = IQC(Psi1,Psi2,M,X,MX_con,[]);
CH=struct();
CH.z=1; CH.w=1; CH.gain="e2p"; CH.gam=-1;
opts = struct();
opts.M1M2_free=true;
info = iqc_analysis(iqc,G,[],CH,opts); % K=[];
disp(info.info);
disp("gam: "+info.gam)

