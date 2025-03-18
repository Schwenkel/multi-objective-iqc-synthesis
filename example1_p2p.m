% Example 1 in [Schwenkel et. al., 2025]
% Revisits Example 14 from [Schwenkel et. al., 2023].
% Therein, a value of gam = 67.81 is achieved as upper bound to the
% p2p-gain using IQCs and a value of 66.93 using pointwise IQCs.
% Define G
A = [0.2 0.01; -.1 -.01];
B = [0.1  0.2 3 2; 
     0.3 -0.2 3 1];
C  = [.2 -.3
      .8 .5
       2 1
       2 3 ];
D   = [ .4 .3   3   1
       -.6 .1   2   7
         1  2   1  -2
        -1  4  -4   3 ];
G = ss(A,B,C,D,-1);
G.InputGroup.p = 1:2;
G.InputGroup.w = 3:4;
G.OutputGroup.q = 1:2;
G.OutputGroup.z = 3:4;
CH = struct();
CH.gain = "p2p";
CH.sigma = 0.6;
CH.w = 1:2;
CH.z = 1:2;
CH.gam=-1;
%% Uncertainty and IQC:
delta1min = -.1;
delta1max = .5;
delta2min = -.3;
delta2max = .6;
ConvHull = [ delta1min delta1min delta1max delta1max;
             delta2min delta2max delta2min delta2max ];
iqc = IQC.polytopic_tv(ConvHull);
%% IQC analysis:
opts.M1M2_free=true;
info1 = iqc_analysis(iqc,G,[],CH,opts); % K=[];
disp(info1.info);
disp("gam with M1, M2 free: "+info1.gam)
%% IQC analysis with M1=sigma*M, M2=(1-sigma)*M, etc. as in (38) in [1]
opts.M1M2_free=false;
info2 = iqc_analysis(iqc,G,[],CH,opts); % K=[];
disp(info2.info);
disp("gam with M1, M2 restricted to (38) with sigma="+CH.sigma+": "+info2.gam)