% Example 2:
% Define G
A = [0.6 0.2; -.1 -.3];
B = [0.2  0.2 3 2 1 0; 
     0.3 -0.2 3 1 2 .2];
C  = [.2 -.3
      .8 .5
       2 1
       2 3
       1 0 ];
D   = [ .4 .3   3   1  0 0
       -.6 .1   2   7  0 .1
         1  2   1  -2  4 0
        -1  4  -4   3  0 0
         0  0  .1  .2  0 0 ];
G = ss(A,B,C,D,-1);
G.InputGroup.p = 1:2;
G.InputGroup.w = 3:4;
G.InputGroup.u = 5:6;
G.OutputGroup.q = 1:2;
G.OutputGroup.z = 3:4;
G.OutputGroup.y = 5;
CH = struct();
CH.gain = "e2e";
CH.sigma = 0.95;
CH.w = 1:2;
CH.z = 1:2;
CH.gam=-1;
CH(2)=CH(1); CH(2).gain = "e2p";
CH(3)=CH(1); CH(3).gain = "p2p";
%% Uncertainty and IQC:
delta1min = -.1;
delta1max = .5;
delta2min = -.3;
delta2max = .6;
rhofilter = -0.25;
nu=4;
iqc1 = IQC.parametric_interval(nu,rhofilter,delta1min,delta1max,1);
iqc2 = IQC.parametric_interval(nu,rhofilter,delta2min,delta2max,1);
iqc = iqc1.combine_with(iqc2);

%% Solve:
K = {[],[],[]};
info = {[],[],[]};
gains = zeros(3);
Ns = [10, 20, 13];
opts = struct();
for i = 1:3
    if CH(i).gain == "p2p"
        opts.rho_start=0.8;
    else
        opts.rho_start=1;
    end
    opts.iterations = Ns(i);
    [Ki,ana_ws,infoi] = iqc_synthesis(G,CH(i),iqc,opts);
    K{i} = Ki;
    info{i} = infoi;
end
%% analyze
for i=1:3
    opts.rho_iter=20;
    ana_wsc = iqc_analysis(iqc,G,K{i},CH,opts);
    for j=1:3
        gains(i,j)=ana_wsc(j).gam;
    end
end
avg_time = zeros(3,1);
gam_lb=zeros(3);
for i = 1:3
    gam_lb(i,:) = lower_bound_exmp2(G,K{i},delta1min,delta1max,delta2min,delta2max);
    for j = 1:length(info{i})
        avg_time(i) = avg_time(i)+info{i}(j).time;
    end
    avg_time(i) = avg_time(i)/length(info{i});
end
%% Display results
disp(" ")
disp("controller | N  | time/N  | resulting e2e         | resulting e2p         | resulting p2p         ")
disp("-----------------------------------------------------------------------------------------------------------")
disp("e2e        | "+ ...
    (length(info{1})-1)+" | "+avg_time(1)+...
    " | "+gam_lb(1,1)+"<=gam<="+gains(1,1)+ ...
    " | "+gam_lb(1,2)+"<=gam<="+gains(1,2)+ ...
    " | "+gam_lb(1,3)+"<=gam<="+gains(1,3))
disp("e2p        | "+ ...
    (length(info{2})-1)+" | "+avg_time(2)+...
    " | "+gam_lb(2,1)+"<=gam<="+gains(2,1)+ ...
    " | "+gam_lb(2,2)+"<=gam<="+gains(2,2)+ ...
    " | "+gam_lb(2,3)+"<=gam<="+gains(2,3))
disp("p2p        | "+...
    (length(info{3})-1)+" | "+avg_time(3)+...
    " | "+gam_lb(3,1)+"<=gam<="+gains(3,1)+ ...
    " | "+gam_lb(3,2)+"<=gam<="+gains(3,2)+ ...
    " | "+gam_lb(3,3)+"<=gam<="+gains(3,3))
%% Comparison to musyn
delta1 = ureal("d1",.3,'Range',[delta1min, delta1max]);
delta2 = ureal("d2",.3,'Range',[delta2min, delta2max]);
Delta = diag([delta1,delta2]);
Gu = lft(Delta,G,2,2);
opts = musynOptions;
opts.MixedMu = 'on';
% musyn uses performance measure gamma that gives the worst-case gain 
% from w to z of the closed loop with 1/gamma*Delta. We want to compare to
% 1*Delta and thus, we rescale the performance channel to find the value
% where musyn provides robust performance = 1
gam_mu_syn = 37.65;
[K{4}, CLperf] = musyn(blkdiag(1/gam_mu_syn*eye(2),1)*Gu,1,2,opts);