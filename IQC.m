classdef IQC < handle
% IQC Class 
% Each iqc contains state space reprentations Psi1 and Psi2, as well as
% a multiplier M and a terminal cost X, which are typically sdpvar and
% restricted to MX_con. If additional sdpvar variables are included in
% MX_con, then they are stacked in the vertical vector H.
%
% The following standard IQCs can be constructed:
%
% iqc = IQC.polytopic_tv(ConvHull)
% iqc = IQC.dynamic(nu,rhofilter,P0)
% iqc = IQC.parametric_interval(nu,rhofilter,alpha,beta,np)
%
properties
    M
    X
    MX_con
    H
    Psi1
    Psi2
end
    
methods

    function self = IQC(Psi1,Psi2,M,X,MX_con,H)
        self.Psi1  = Psi1;
        self.Psi2  = Psi2;
        self.M   = M;
        self.X   = X;
        self.MX_con = MX_con;
        self.H   = allvariables(H);
    end
    
    function Psi_ss = Psi(self)
        Psi_ss = [self.Psi1, self.Psi2];
    end

    function copy = copy(self)
        copy = IQC(self.Psi1,self.Psi2,self.M,self.X,self.MX_con,self.H);
    end

    function n = n1(self)
        n = length(self.Psi1.A);
    end

    function n = n2(self)
        n = length(self.Psi2.A);
    end
    
    function [M2, X2, MX2_con, H2] = get_MX2(self)
        % makes a copy of the sdp variables M, X, and H as well as the
        % constraints MX_con.
        h1 = allvariables(self.M);
        h2 = allvariables(self.X);
        h3 = allvariables(self.H);
        hold = [h1;h2;h3];
        hnew = sdpvar(length(hold),1);
        M2 = replace(self.M,hold,hnew);
        X2 = replace(self.X,hold,hnew);
        H2 = replace(self.H,hold,hnew);
        MX2_con = replace(self.MX_con{1},hold,hnew);
        MX2_con = {MX2_con, self.MX_con{2}};
    end
    
    function comb = combine_with(self, other)
        % stack IQCs 
        cPsi1   = blkdiag(self.Psi1, other.Psi1);
        cPsi2   = blkdiag(self.Psi2, other.Psi2);
        cM      = blkdiag(self.M,    other.M);
        sX11    = self.X(1:self.n1, 1:self.n1);
        sX12    = self.X(1:self.n1, self.n1+1:end);
        sX22    = self.X(self.n1+1:end, self.n1+1:end);
        oX11    = other.X(1:other.n1, 1:other.n1);
        oX12    = other.X(1:other.n1, other.n1+1:end);
        oX22    = other.X(other.n1+1:end, other.n1+1:end);
        cX      = [blkdiag(sX11,oX11)   blkdiag(sX12,oX12)
                   blkdiag(sX12,oX12)'  blkdiag(sX22,oX22)];
        cMX_con = {[self.MX_con{1} other.MX_con{1}], ...
                  [self.MX_con{2}; other.MX_con{2}]};
        cH      = [self.H; other.H];
        comb = IQC(cPsi1,cPsi2,cM,cX,cMX_con,cH);
    end
end

methods (Static)
    
    function iqc_p_i = parametric_interval(nu,rhofilter,dmin,dmax,np)
        % Returns IQC for scalar parametric uncertainties p = delta*q with
        % delta in the interval [dmin, dmax]
        % np = dimension of p
        % nu = dimension of IQC filter
        % rhofilter = location of poles of IQC filter
        % This IQC is the discrete-time version of the IQC with terminal 
        % cost introduced in [Scherer, Veenman, 2018].
        %
        % auxiliary variable
        H = sdpvar(2*nu);
        % terminal cost
        K = sdpvar(nu, nu, 'full');
        X = [zeros(nu) K'; K zeros(nu)];
        % multiplier
        nP = np*(nu+1);
        P = sdpvar(nP, nP, 'full');
        M = [zeros(nP) P'; P zeros(nP)];
        % The Filter
        phi = IQC.get_phi(nu,rhofilter);     % scalar
        % phi = kron(eye(np),phi); % kron does not work with ss, hency by 
        % hand:
        Ap = phi.A;                  % p dimensional
        Bp = repmat(phi.B, 1,np);    % p dimensional
        Cp = repmat(phi.C,np,1);     % p dimensional    
        Dp = kron(eye(np),phi.D);    % p dimensional
        Psi1 = ss(Ap,Bp,[dmax*Cp;-dmin*Cp],[dmax*Dp;-dmin*Dp],-1);
        Psi2 = ss(Ap,Bp,[-Cp;Cp],[-Dp;Dp],-1);
        Psi = [Psi1 Psi2];
        % We transformed ssrep of Psi from Scherer, Veenman (2018) by T to 
        % have the explicit minimal representation of Psi1 and Psi2. 
        % Thus we must transform the terminal cost by T as well:
        T = [dmax*eye(nu) -eye(nu); -dmin*eye(nu) eye(nu)];
        H = T'*H*T; % makes no difference in theory as H=H' has no further 
                    % structure but seems to be numerically more stable.
        X = T'*X*T;
        % The LMI constraints:
        LMI1 = H-X<=0; eps1=0;
        outer1 = [eye(2*nu)  zeros(2*nu, np); 
                  Psi.A      Psi.B(:,1:np); 
                  Psi.C      Psi.D(:,1:np)];
        inner1 = blkdiag(-H, H, M);
        LM2 = outer1'*inner1*outer1;
        eps2 = 1e-6;
        LMI2 = LM2 >= eps2*eye(length(LM2));
        MX_con = {[LMI1 LMI2], [eps1; eps2]};
        iqc_p_i = IQC(Psi1,Psi2,M,X,MX_con,H);
    end

    function iqc_d = dynamic(nu,rhofilter,P0)
        % Returns IQC for SISO dynamic uncertainties p = delta*q where
        % delta satisfies the FDI [1; delta]'*P0*[1;delta] > 0.
        % nu = dimension of IQC filter
        % rhofilter = location of poles of IQC filter
        % This IQC is the discrete-time version of the IQC with terminal 
        % cost introduced in [Scherer, 2022].
        %
        % terminal cost
        K = sdpvar(nu);
        X = kron(P0,K);
        m = sdpvar(nu+1);
        M = kron(P0,m);
        % The Filter
        phi = IQC.get_phi(nu,rhofilter);
        Psi1 = ss(phi.A,phi.B,[phi.C; 0*phi.C],[phi.D; 0*phi.D],-1);
        Psi2 = ss(phi.A,phi.B,[0*phi.C; phi.C],[0*phi.D; phi.D],-1);
        % The LMI constraints:
        outer1 = [eye(nu)  zeros(nu,1); 
                  phi.A    phi.B; 
                  phi.C    phi.D ];
        inner1 = blkdiag(-K, K, m);
        MX_con = outer1'*inner1*outer1 >= 0;
        MX_con = {MX_con, 0};
        H=[]; % no auxiliary variables
        iqc_d = IQC(Psi1,Psi2,M,X,MX_con,H);
    end
        
    function iqc_p_tv = polytopic_tv(ConvHull)
        % time-varying polytopic uncertainty diag(delta) with delta in the
        % convex hull of the column vectors of ConvHull.
        % This IQC can be found in [Schwenkel et. al., 2023] and is the 
        % discrete time version of "Class 4" in [Veenman et. al., 2016].
        %
        [n,d] = size(ConvHull);
        A = [];
        Dq = [eye(n); zeros(n)];
        Dp = [zeros(n); eye(n)];
        B = zeros(0, n);
        C = zeros(2*n, 0);
        Psi1 = ss(A,B,C,Dq);
        Psi2 = ss(A,B,C,Dp);
        M = sdpvar(2*n);
        X=[];
        lambda_con = Dp'*M*Dp <= 0;
        margins = 0;
        for i=1:d
            ID = [eye(n); diag(ConvHull(:,i))];
            lambda_con = [lambda_con, -ID'*M*ID <=0];
            margins = [margins;0];
        end
        iqc_p_tv = IQC(Psi1, Psi2, M, X, {lambda_con,margins}, []);
    end

    function phi = get_phi(nu,rho)
        % returns ss of phi = [ 1; 1/(z-rho); ... ; 1/(z-rho)^nu ] 
        % from (12a) and (21) from [Veenman et. al., 2016].
        %
        A = rho*eye(nu)+diag(ones(nu-1,1),-1);
        B = [1; zeros(nu-1,1)];
        C = [zeros(1, nu); eye(nu)];
        D = [1; zeros(nu,1)];
        phi=ss(A,B,C,D,-1);
    end
        
    function iqc_combined = combine_multiple(varargin)
        iqc_combined = varargin{1};
        for i=2:nargin
            iqc_combined = iqc_combined.combine(varargin{i});
        end
    end
end
    
end