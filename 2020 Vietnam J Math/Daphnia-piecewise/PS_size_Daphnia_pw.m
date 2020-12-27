function out = PS_size_Daphnia_pw
% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Breda, Diekmann, Gyllenberg, Vermiglio (2020), Numerical
% bifurcation analysis of physiologically structured population models via
% pseudospectral approximation, Vietnam J Math
%
%% PS_size_Daphnia_pw.m
% Matcont system definition file of the PseudoSpectral Discretization of
% the Daphnia model with rates as in Breda et al, SISC 2015
% Discretization of size, piecewise approach

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= @userf; %[];
out{11}= []; %@userf2;
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,K,mu,aux,M) 
% Parameters MUST be listed separately

    xb=0.8; xA=2.5; xm=6; gamma=0.15; xi=7; nu=1.8; rm=0.1; 
    a1=0.5;%a1=2;%

c1 = state(1:M);
c2 = state(M+1:2*M);

%% Discretization of the interval [xb,xA] and [xA,xm] (piecewise)

% juvenile discretization
[QuadWeights_j,Nodes_j,DD_j,BaryWeights_j]=cheb(M,xb,xA);
DM_j = DD_j(2:end,2:end);

% adult discretization
[QuadWeights_a,Nodes_a,DD_a,BaryWeights_a]=cheb(M,xA,xm);
DM_a = DD_a(2:end,2:end);
    
%% SYSTEM DEFINITION *** to be completed by the user ***
% DAPHNIA MODEL

% Parameters and functions

    % resource intrinsic rate of change
    free_dyn=@(S) a1*S*(1-S/K); % logistic
    %free_dyn= @(S) alpha*(K-S); % chemostat
    
    % Holling type II functional response
    func_resp = @(S) xi*S./(1+xi*S);
    BBeta= @(X,S) rm*func_resp(S).*(X.^2).*(X>=xA);
    GGamma= @(X,S) nu*func_resp(S).*(X.^2);
   
% growth rate
g= @(X,S) max(0,gamma*(xm*xi*S./(1+xi*S)-X)); 
% g0=0.1;
% g= @(X,S) max(g0,gamma*(xm*xi*S./(1+xi*S)-X)); 

% Ell_prime = zeros(M+1,M);
% for jj=1:M
%     Ell_prime(:,jj) = interpoly(Nodes_Adults,Nodes,DD(:,jj+1),BaryWeights);
% end
% 
% tildeBeta = (QuadWeights_Adults*(BBeta(Nodes_Adults,state(end)).*Ell_prime));
% tildeGamma =(QuadWeights*(GGamma(Nodes,state(end)).*DD(:,2:end)));

Beta_interp = zeros(M+1,M+1); 
for jj=1:M+1
    Beta_interp(:,jj) = BBeta(Nodes_a,state(end)).*DD_a(:,jj);
end

Gamma_interp_j = zeros(M+1,M+1);
Gamma_interp_a = zeros(M+1,M+1);

for jj=1:M+1
    Gamma_interp_j(:,jj) = GGamma(Nodes_j,state(end)).*DD_j(:,jj);
    Gamma_interp_a(:,jj) = GGamma(Nodes_a,state(end)).*DD_a(:,jj);
end

tildeBeta = (QuadWeights_a*Beta_interp);
tildeGamma_j =(QuadWeights_j*Gamma_interp_j);
tildeGamma_a =(QuadWeights_a*Gamma_interp_a);

GG = diag(g([Nodes_j(2:end);Nodes_a(2:end)],state(end)));
Sigma = mu*eye(2*M);

%% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION

DM = [ DM_j, zeros(M,M);
       zeros(M,M-1), DD_a(2:end,:)];
%blkdiag(DM_j,DM_a);
    
    dydt= [
		- GG*DM*[c1;c2]+(tildeBeta*[c1(end);c2])*ones(2*M,1)-Sigma*[c1;c2] + (1-aux);
		free_dyn(state(end))-tildeGamma_j*[0;c1]-tildeGamma_a*[c1(end);c2]+(1-aux)
        ];
    
end
 
 
% --------------------------------------------------------------------------
function Weq=init(M,xeq,yeq)
% INPUT: M discretization index 
%        xeq,yeq scalar values of equilibria
% OUTPUT Weq is the initial vector for init_EP_EP
    
    xb=0.8; xA=2.5; xm=6; 
% juvenile discretization
    [~,Nodes_j,~,~]=cheb(M,xb,xA);

% adult discretization
    [~,Nodes_a,~,~]=cheb(M,xA,xm);

Weq= [xeq*(Nodes_j(2:end)-Nodes_j(1)); xeq*(Nodes_a(2:end)-Nodes_j(1)); yeq];
 
end

function y = userf(time,state,K,mu,aux,M) 
    y = (K-1)*(K-1.5)*(K-2);
end

%% AUXILIARY FUNCTIONS (Chebyshev discretization & barycentric interpolation)

function [w,x,D,q]=cheb(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% D - differentiation matrix
% q - row vector of the barycentric weights
% see Trefethen

if N==0
    x=1;
    D=0;
    return
end
p=pi*(0:N)'/N;
x=((a-b)*cos(p)+b+a)/2; % ordered from a to b
% x=((b-a)*cos(p)+b+a)/2; % ordered from b to a

c=[2;ones(N-1,1);2].*(-1).^(0:N)';
X=repmat(x,1,N+1);%X=x(:,ones(1,N+1)); 
dX=X-X';
D=(c*(1./c)')./(dX+(eye(N+1)));
D=D-diag(sum(D,2)); %D=D-diag(sum(D'));

% Quadrature weights
w=zeros(1,N+1);
ii=2:N;
v=ones(N-1,1);
if mod(N,2)==0
    w(1)=1/(N^2-1);
    w(N+1)=w(1);
    for k=1:N/2-1
        v=v-2*cos(2*k*p(ii))/(4*k^2-1);
    end
    v=v-cos(N*p(ii))/(N^2-1);
else
    w(1)=1/N^2;
    w(N+1)=w(1);
    for k=1:(N-1)/2
        v=v-2*cos(2*k*p(ii))/(4*k^2-1);
    end
end
w(ii)=2*v/N;
w=w*abs(b-a)/2;

% Barycentric weights
q=1./prod(dX'+eye(N+1)); %q=1./prod(dX'+eye(N+1)); % row vector of the barycentric weights
end
