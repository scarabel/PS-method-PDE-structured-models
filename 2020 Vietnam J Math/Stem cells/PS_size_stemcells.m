function out = PS_size_stemcells
% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Breda, Diekmann, Gyllenberg, Vermiglio (2020), Numerical
% bifurcation analysis of physiologically structured population models via
% pseudospectral approximation, Vietnam J Math
%
%% PS_size_stemcells.m
% Matcont system definition file of the PseudoSpectral Discretization of
% the PDE for Getto et al JMB 2019 and Doumic et al 2011

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= @userf;
out{11}= [];
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,a,p,ka,kp,mu,muw,x1,x2,aux,M) 
% Parameters MUST be listed separately

%% Discretization of the size interval [x1,x2]

    [QuadWeights,Nodes,DD,BaryWeights]=cheb(M,x1,x2);
	DM = DD(2:end,2:end);


%% SYSTEM DEFINITION *** to be completed by the user ***
% STEM CELL MODEL
% The state vector is (c,w,v) where c is progenitors, w stem cells, v
% adults

% Parameters and functions

    dw = @(v) p./(1+kp*v);
    sw = @(v) a./(1+ka*v);
    qq = @(v) (2*sw(v)-1).*dw(v)-muw;
    gamma = @(v) 2*(1-sw(v)).*dw(v); 
    
    gg = @(mat,v) 2*p*(1-a./(1+v));

    % Discretization 
    
    WW = state(end-1); VV = state(end); 
    GG = diag(gg(Nodes(2:end),VV));
    
%% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION
 
    select_last(M) = 1; % used to select the last component of a column vector

    dydt= [
		- GG*DM*state(1:M) + gamma(VV)*WW*ones(M,1) + (1-aux);
		qq(VV)*WW + (1-aux);
        gg(x2,VV)*select_last*(DM*state(1:M)) - mu*VV
        ];
    
end
 
 
% --------------------------------------------------------------------------
function Eq=init(M,meq,Weq,Veq,x1,x2)
% INPUT: M discretization index 
%        meq,Weq,Veq scalar values of equilibria
% OUTPUT Weq is the initial vector for init_EP_EP
    
    [~,Nodes,~,~] = cheb(M,x1,x2);
    Eq = [meq*(Nodes(2:end)-Nodes(1));Weq;Veq];
 
end

function y = userf(time,state,a,p,ka,kp,mu,muw,x1,x2,aux,M) 
    y = (p-3.9); %(p-2);%
end

%% AUXILIARY FUNCTIONS (Chebyshev discretization)

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
