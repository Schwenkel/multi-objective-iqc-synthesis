function gam_lb = lower_bound_exmp2(G,K,delta1min,delta1max,delta2min,delta2max)
%% find worst case p2p and e2e input
N = 2;
w = G.InputGroup.w;
nw = length(w);
[theta1, theta2] = meshgrid(linspace(0,1,N),linspace(0,1,N));
%theta1 = 0*theta1;
%theta2 = 0*theta2;
gamma_p2p = 0*theta1;
gamma_e2e = 0*theta1;
GK = lft(G,K);
nx = length(GK.A);
for i=1:N
    for j=1:N
        delta = [theta1(i,j)*delta1min+(1-theta1(i,j))*delta1max
                 theta2(i,j)*delta2min+(1-theta2(i,j))*delta2max];
        Deltass = ss([], zeros(0,2), zeros(2,0), [delta(1) 0;0 delta(2)],-1);
        CL = lft(Deltass, GK, 2, 2);
        gamma_e2e(i,j)=hinfnorm(CL); 
        T = 30;
        w = sdpvar(nw,T);
        assign(w,ones(nw,T)/sqrt(nw))
        LMI = [];
        xk = zeros(nx,1);
        for k=1:T-1
            LMI=[LMI, [eye(nw) w(:,k); w(:,k)' 1] >=0];
            xk = CL.A*xk+CL.B*w(:,k);
        end
        LMI=[LMI, [eye(nw) w(:,T); w(:,T)' 1] >=0];
        zT = CL.C*xk+CL.D*w(:,T);
        gam = sdpvar(1);
        options = sdpsettings('solver','mosek','verbose',0);
        for k=1:4
            zT_guess = value(zT);
            LMI2 = 2*zT'*zT_guess-zT_guess'*zT_guess-gam >=0;
            sol = optimize([LMI, LMI2], -gam, options);
            if sol.problem
                error("LMI not solved")
            end
        end
        gamma_p2p(i,j)= sqrt(value(gam));
    end
end

%% find worst case e2p input
N = 1;
[theta1, theta2] = meshgrid(linspace(0,1,N),linspace(0,1,N));
theta1 = 0*theta1;
theta2 = 0*theta2;
gamma_e2p = 0*theta1;
for i=1:N
    for j=1:N
        delta = [theta1(i,j)*delta1min+(1-theta1(i,j))*delta1max
                 theta2(i,j)*delta2min+(1-theta2(i,j))*delta2max];
        Deltass = ss([], zeros(0,2), zeros(2,0), [delta(1) 0;0 delta(2)],-1);
        CL = lft(Deltass, GK, 2, 2);
        T = 20;
        w = sdpvar(nw,T);
        assign(w,ones(nw,T)/sqrt(nw))
        LMI = [eye(nw*T) w(:); w(:)' 1] >=0;
        xk = zeros(nx,1);
        for k=1:T-1
            xk = CL.A*xk+CL.B*w(:,k);
        end
        zT = CL.C*xk+CL.D*w(:,T);
        gam = sdpvar(1);
        options = sdpsettings('solver','mosek','verbose',0);
        for k=1:4
            zT_guess = value(zT);
            LMI2 = 2*zT'*zT_guess-zT_guess'*zT_guess-gam >=0;
            sol = optimize([LMI, LMI2], -gam, options);
            if sol.problem
                error("LMI not solved")
            end
        end
        gamma_e2p(i,j)= sqrt(value(gam));
    end
end
gam_lb(1) = max(gamma_e2e(:));
gam_lb(2) = max(gamma_e2p(:));
gam_lb(3) = max(gamma_p2p(:));