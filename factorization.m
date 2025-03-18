function [Psih, Xh, Vh, Z] = factorization(Psi1,Psi2,M,X)
% returns a factorization 
% Psih'*[I 0; 0 -I]*Psih = [Psi1 Psi2]'*M*[Psi1 Psi2],
% the terminal cost Xh in the transformed filter state, as well as the 
% filter state transformation matrix Vh, and the certificate Z.
% The factorization has the properties of Theorem 4 in 
% [Schwenkel et. al. 2025] and is an implementation of the constructive
% proof of this Theorem. 
nq = size(Psi1.B,2);
np = size(Psi2.B,2);

% 1. Constructing Psih_11
Q = Psi1.C'*M*Psi1.C;
R = Psi1.D'*M*Psi1.D;
S = Psi1.C'*M*Psi1.D;
[Zu, ~, n0] = dare_anti(Psi1.A,Psi1.B,Q,R,S);
D11hatu = chol(Psi1.B'*Zu*Psi1.B+R);
C11hatu = D11hatu'\(Psi1.B'*Zu*Psi1.A+S');
Psi11u = ss(Psi1.A,Psi1.B,C11hatu,D11hatu,-1);
Psi11 = tf('z')^(-n0)*Psi11u;
Psi11Psi1 = minreal([Psi11; Psi1],[],false);
Psi11Psi1 = balreal(Psi11Psi1);
A1hat = Psi11Psi1.A;              B1hat = Psi11Psi1.B; 
C11hat = Psi11Psi1.C(1:nq,:);     D11hat = Psi11Psi1.D(1:nq,:);   
C1hat = Psi11Psi1.C(nq+1:end,:);

% 2. Constructing Psih_12
B1inv = B1hat/D11hat;
D1inv = Psi1.D/D11hat;
Psi1Psi11inv = ss(A1hat-B1inv*C11hat,B1inv,C1hat-D1inv*C11hat,D1inv,-1);
[~, Psi11Psi11invMPsi2] = isproper(Psi1Psi11inv'*M*Psi2);
Psi12Psi2 = minreal([Psi11Psi11invMPsi2; Psi2],[],false);
Psi12Psi2 = balreal(Psi12Psi2);
A2hat = Psi12Psi2.A;              B2hat = Psi12Psi2.B; 
C12hat = Psi12Psi2.C(1:nq,:);     D12hat = Psi12Psi2.D(1:nq,:);   
C2hat = Psi12Psi2.C(nq+1:end,:);

% 3. Constructing Psih_22
Q2 = Psi12Psi2.C'*blkdiag(eye(nq),-M)*Psi12Psi2.C;
R2 = Psi12Psi2.D'*blkdiag(eye(nq),-M)*Psi12Psi2.D;
S2 = Psi12Psi2.C'*blkdiag(eye(nq),-M)*Psi12Psi2.D;
Q2 = (Q2+Q2')/2;
R2 = (R2+R2')/2;
Zs = idare(A2hat,B2hat,Q2,R2,S2);
D22hat = chol(B2hat'*Zs*B2hat+R2);
C22hat = D22hat'\(B2hat'*Zs*A2hat+S2');

% 4. Defining Psih
Ahat = blkdiag(A1hat,A2hat);
Bhat = blkdiag(B1hat,B2hat);
Chat = [C11hat                   C12hat;
        zeros(np,length(A1hat))  C22hat];
Dhat = [D11hat        D12hat;
        zeros(np,nq)  D22hat];
Psih.A1 = A1hat;
Psih.A2 = A2hat;
Psih.A = Ahat;
Psih.B1 = B1hat;
Psih.B2 = B2hat;
Psih.B = Bhat;
Psih.C1 = C11hat;
Psih.C2 = C22hat;
Psih.C3 = C12hat;
Psih.C = Chat;
Psih.D1 = D11hat;
Psih.D2 = D22hat;
Psih.D3 = D12hat;
Psih.D = Dhat;
Psih.Chat = [C1hat C2hat];
Psih.Dhat = [Psi1.D Psi2.D];

% 5. Computing Z, Xh, and Vh
Mhat = blkdiag(eye(nq),-eye(np));
Q = [Chat; C1hat C2hat]'*blkdiag(Mhat,-M)*[Chat; C1hat C2hat];
Q = (Q+Q')/2; % ensure symmetry
Z = dlyap(Ahat',Q);
[Xh, Vh] = compute_Xhat(Psi1, Psi2, Psih, Z, X);

end

function [Xh, Vh] = compute_Xhat(Psi1, Psi2, Psih, Z, X)
    % find Vhat such that Vhat*Ahat=Apsi*Vhat, Vhat*Bhat=Bpsi, and
    % Chat=Cpsi*Vhat by solving system of linear equations.
    npsih = length(Psih.A1)+length(Psih.A2);
    Apsi = blkdiag(Psi1.A, Psi2.A);
    Bpsi = blkdiag(Psi1.B, Psi2.B);
    Cpsi = [Psi1.C Psi2.C];
    npsi = length(Apsi);
    Vh = [ kron(eye(npsih),Apsi)-kron(Psih.A',eye(npsi));
             kron(eye(npsih),Cpsi);
             kron(Psih.B',eye(npsi))                        ] \ ...
           [ zeros(npsi*npsih,1); Psih.Chat(:); Bpsi(:)       ];
    Vh = reshape(Vh,[npsi, npsih]);
    % Alternative way to compute Vhat
    % [Ah,~,Ch,T1] = obsvf(Psih.A, Psih.B, Psih.Chat);
    % T2 = obsv(Ah(nW+1:end,nW+1:end), Ch(:,nW+1:end));
    % T3 = obsv(Psi.A, Psi.C);
    % Vhat = T3\T2*T1(nW+1:end,:);
    Xh = Vh'*X*Vh+Z;
    Xh = (Xh+Xh')/2;
end

function [X,K,n_zero] = dare_anti(A,B,Q,R,S)
    % computes the solution X of the dare
    % A'*X*A-X-(A'*X*B+S)*((B'*X*B+R)\(A'*X*B+S)')+Q = 0
    % such that eig(A+B*K) has n_zero eigenvalues at 0 and all others between
    % 1 and Inf, i.e., 1 < |lambda_i| < Inf
    %
    % Based on the procedure from [Ionescu, Weiss, 1992].
    %
    nx = length(A);
    nu = size(B,2);
    % scaling for numerical improvements:
    scale = max(max(abs([Q(:); R(:); S(:)])))/100;
    Q = Q/scale;
    R = R/scale;
    S = S/scale;
    % define extended simplectic pencil (ESP) as in (11) in [1]
    M = [eye(nx)      zeros(nx)  zeros(nx,nu); 
         zeros(nx)    -A'        zeros(nx,nu); 
         zeros(nu,nx) -B'        zeros(nu,nu)];
    N = [A   zeros(nx)     B; 
         Q   -eye(nx)      S; 
         S'  zeros(nu,nx)  R];
    % separate the nu first order infinite eigenvalues (M has always rank 2*nx) 
    [U,MS,Z1] = svd(M);
    U = U';
    NS = U*N*Z1;
    i1 = 1:2*nx;
    i2 = 2*nx+1:2*nx+nu;
    % Ensure the nu eigenvalues at Inf are really Inf:
    MS(i2,i2) = 0;
    % Make NS(i2,i1)=0 and NS(i2,i2) an upper triangular matrix 
    [Z2, NS2] = qr(NS(i2,:)');      % Now NS2' is lower triangular
    NS2 = NS2';                     % Now NS(i2,:) = NS2*V2'
    NS2 = flip(flip(NS2,1),2);      % [ * 0 0; * * 0]  -->  [ 0 * *; 0 0 *]
    Z2 = flip(Z2,2);                % Now U*N*Z1*Z2 = [ * * *; 0 0 *; 0 * *]
    % U is not needed, but for completeness:
    % U(i2,:) = flip(U(i2,:),1);    % Now U*N*Z1*Z2 = [ * * *; 0 * *; 0 0 *]
    % Define Z=Z1*Z2 and NS = U*N*Z; MS = U*N*Z, i.e.
    Z = Z1*Z2;
    NS(i1,:) = NS(i1,:)*Z2;
    NS(i2,:) = NS2;
    MS(i1,:) = MS(i1,:)*Z2;
    % MS(i2,:) = flip(MS(i2,:),2);  % just flips around zeros
    % Now we have
    % NS = U*N*Z = [ NS11 NS12 ]
    %              [  0   NS22 ]
    % MS = U*M*Z = [ MS11 MS12 ]
    %              [  0    0   ]
    % with NS22 non-singular and upper triangular
    % Next, we bring lambda*MS11-NS11 to its generalized Schur decomposition:
    [NN,MM,QQ,ZZ] = qz(NS(i1,i1),MS(i1,i1),'real');
    lambda = abs(ordeig(NN,MM));               % all eigenvalues
    clusters = 2*ones(2*nx,1);                 % stable cluster (w/o zeros)
    i_inf = find(lambda == Inf);               % index of infinite eigenvalues
    [~,i_zero] = mink(lambda,length(i_inf));   % index of corresponding zeros
    clusters(i_zero) = 4;                      % zero eigenvalue cluster
    clusters(i_inf) = 1;                       % infinite eigenvalue cluster
    clusters(lambda ~= Inf & lambda > 1) = 3;  % anti-stable cluster
    if sum(clusters>=3) ~= nx
        % this typically occurs if Assumption 1 is not satisfied.
        error("something went wrong with separating the eigenvalues")
    end
    [~,~,~,ZZS] = ordqz(NN,MM,QQ,ZZ,clusters);
    Z(:,i1) = Z(:,i1)*ZZS;
    % U is not needed, but for completeness:
    % U(i1,:) = QQS*U(i1,:); % where [~,~,QQS,~] = ordqz(NN,MM,QQ,ZZ,clusters);
    % Now U*M*Z is triu and U*N*Z is blktriu with blocks of dim 2 for complex eigenvalues;
    i1a = 1:nx;
    i1b = nx+1:2*nx;
    V1 = Z(i1a,i1a);
    V2 = Z(i1b,i1a);
    V3 = Z(i2,i1a);
    X = V2/V1;
    X = scale*X;
    K = V3/V1;
    n_zero=length(i_zero);
end

