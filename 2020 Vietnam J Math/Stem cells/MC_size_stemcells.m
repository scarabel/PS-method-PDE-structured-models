% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Breda, Diekmann, Gyllenberg, Vermiglio (2020), Numerical
% bifurcation analysis of physiologically structured population models via
% pseudospectral approximation, Vietnam J Math
%
%% MC_size_stemcells
% command line instructions for Matcont continuation of the system defined
% in PS_size_stemcells

%profile on
clear; close all

M=10; 
MM= M+2;

% initial value of continuation parameter
a=0.9;
p=1; % p=21, mu=100, epsilon=1; muw=1;
ka=1;
kp=0;
mu=5;
muw=1; 
x1=1;
x2=2;
aux=1;


ap1=2; % index of the bifurcation parameter in the vector 'par'
ap2=5;

% case gg = @(mat,v) 2*p*(1-a./(1+v));
% (s)s
ka=1; kp=0; mu=1.75; p=1.5; ap1=2; ap2=5; figure(2); axis([0 2.5 1.9 5.1]); figname=[num2str(M),'_ss_hopf']; matname=[num2str(M),'_ss']; % (s)s
%ka=0; kp=1; mu=8; p=1.4; ap1=2; ap2=5; figure(2); axis([3 8 1 2.5]); figname=[num2str(M),'_sp_hopf']; matname=[num2str(M),'_sp']; % (s)p



par=[a,p,ka,kp,mu,muw,x1,x2,aux,M]';


% initial equilibrium point

% tolerance
TOL=1e-3;
TestTOL=TOL;
VarTOL=TOL;

%% Continuation

handles=feval(@PS_size_stemcells); 
opt=contset; 
global cds

%% Initilize equilibrium vector

rhs = @(t,y) feval(handles{2},t,y,a,p,ka,kp,mu,muw,x1,x2,aux,M);
[TOUT,YOUT] = ode45(rhs,[0 500],[ones(M,1);0.1;0.1]);
figure(5)
plot(TOUT,YOUT(:,end))

Eq = YOUT(end,:)';

%% Equilibrium continuation from initial point [xeq;yeq]
% xeq = equilibrium of RE
% yeq = equilibrium of DDE
% par = vectors of parameters

disp('Starting equilibrium continuation from initial point');
par0=par;

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'MaxStepSize',0.1); 
opt=contset(opt,'Backward',0);

%     % detection of the periodic orbits
%     UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P';
%     opt=contset(opt,'Userfunctions',1);
%     opt=contset(opt,'UserfunctionsInfo',UserInfo);
    
%Eq = feval(handles{1},M,Ueq,Weq,Veq); % initializes equilibrium vector

[x0,v0]=init_EP_EP(@PS_size_stemcells,Eq,par0,ap1);

tic
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt);
time1=toc;

xe(end,end)


%% Detection of singular points
% H, Hopf point

for ii=size(se):-1:1
    if strcmp(se(ii).label,'H ')==1
        H_index=se(ii).index;
        break;
    end
end
par(ap1)=xe(end,H_index);
H=xe(1:MM,H_index);

xeH=xe; veH=ve; seH=se; heH=he; feH=fe;
parH=par;


% P1, P2 corresponding to p=3.9
for ii=1:length(se)
    if (strcmp(se(ii).label,'P')==1 && abs((xe(end,se(ii).index)-3.9))<0.01)
        ind_eq39=se(ii).index;
    end
end

%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

% set options
TOL= 1e-3;
opt=contset(opt,'FunTolerance',TOL); 
opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);

opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'Eigenvalues',0);
opt=contset(opt,'Backward',0);
opt=contset(opt,'MaxStepsize',0.1);


[x0,v0]=init_H_H(@PS_size_stemcells,H,parH,[ap2 ap1]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)
[xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)

jj=0;
while (xh(MM+1,end)>0 && xh(MM+1,end)<8 && jj<10)
      [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
      jj=jj+1;
end

% Plot
figure(2)
cpl(xh,vh,sh,[MM+1 MM+2]); hold on
%xlabel('mu'); ylabel('K');

%%
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',0);
opt=contset(opt,'Backward',1);
opt=contset(opt,'MaxStepsize',0.1);

%opt=contset(opt,'MaxNumPoints',15);
%opt=contset(opt,'Adapt',0);

[x0,v0]=init_H_H(@PS_size_stemcells,H,parH,[ap2 ap1]);
[xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(MM+1,end)
[xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
jj=0;
while (xhb(MM+1,end)>0 && xhb(MM+1,end)>0 && jj<10)
      [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
      jj=jj+1;
end

cpl(xhb,vhb,shb,[MM+1 MM+2]); hold on
title(['Hopf curve, M=',num2str(M)]);
xlabel('mu')
ylabel('p')

% savefig(figname);


%% Limit cycle continuation from H
% H = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from H');

% set options
opt=contset(opt,'Backward',0);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxStepsize',0.1);
% opt=contset(opt,'InitStepsize',1e-4);
% opt=contset(opt,'MinStepsize',1e-6);
opt=contset(opt,'MaxNumPoints',50);
%opt=contset(opt,'Adapt',0);

%TOL=1e-3;
%opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
%opt=contset(opt,'TestTolerance',1e-2);

%opt=contset(opt,'Increment',1e-6); %default 1e-5
%opt=contset(opt,'MaxCorrIters',100); %default 1e-5
%opt=contset(opt,'MaxNewtonIters',10); %default 1e-5

    % detection of the periodic orbits
    UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P';
    opt=contset(opt,'Userfunctions',1);
    opt=contset(opt,'UserfunctionsInfo',UserInfo);
    
ntst=20; % number of interval
ncol=2; % degree of polynomial

[x0,v0]=init_H_LC(@PS_size_stemcells,H,parH,ap1,0.01,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
[xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
[xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)

%%
% jj=0;
% while (max(xlc(end,:))<4 && jj<5)
%     [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
%     jj=jj+1;
% end

% Plot max and min periodic solutions
    % re-interpolation for a smoother plot
    mesh_refined=linspace(0,1,100);
    Per_Solutions_m = zeros(length(mesh_refined),size(xlc,2));
    Per_Solutions_w = zeros(length(mesh_refined),size(xlc,2));
    Per_Solutions_v = zeros(length(mesh_refined),size(xlc,2));
    for ind_persol=1:size(xlc,2)
%         size(flc(1:ntst+1,ind_persol))
%         size(xlc(1:MM:((ntst*ncol+1)*MM),ind_persol))
%         size(mesh_refined)
        Per_Solutions_m(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(M:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        Per_Solutions_w(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(M+1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        Per_Solutions_v(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(M+2:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
    end
    
upperbound_m=max(Per_Solutions_m);
lowerbound_m=min(Per_Solutions_m);
upperbound_w=max(Per_Solutions_w);
lowerbound_w=min(Per_Solutions_w);
upperbound_v=max(Per_Solutions_v);
lowerbound_v=min(Per_Solutions_v);

index=find(min(abs(flc(ntst+2:end,:)-1))<0.1,1,'last') %size(xlc,2);

figure
subplot(3,1,1)
cpl(xe,ve,se,[MM+1 M]); hold on
plot(xlc(end,1:index),upperbound_m(1:index),'g',xlc(end,1:index),lowerbound_m(1:index),'g');
xlabel('p','interpreter','latex');
title('bifurcation of m(t,xbar)');
%axis([0 2 0 0.05])
%ylabel('max/min','interpreter','latex')

subplot(3,1,2)
cpl(xe,ve,se,[MM+1 M+1]); hold on
plot(xlc(end,1:index),upperbound_w(1:index),'g',xlc(end,1:index),lowerbound_w(1:index),'g');
xlabel('p','interpreter','latex');
title('bifurcation of w');
%axis([0 2 0 3])

subplot(3,1,3)
cpl(xe,ve,se,[MM+1 M+2]); hold on
plot(xlc(end,1:index),upperbound_v(1:index),'g',xlc(end,1:index),lowerbound_v(1:index),'g');
xlabel('p','interpreter','latex');
title('bifurcation of v');

%figname=[num2str(M),'_bif_stemcells'];
%savefig(figname);

%% Plot profile periodic orbit in real period

% orbit (unstable and stable) in p=3.9

ind_persol1=[];
for ii=1:length(slc)
    if (strcmp(slc(ii).label,'P')==1 && isempty(ind_persol1))
        ind_persol1=slc(ii).index;
    elseif (strcmp(slc(ii).label,'P')==1 && ~isempty(ind_persol1))
        ind_persol2=slc(ii).index;
    end
end

% ind_persol1=find(xlc(end,:)>1,1,'first');
% ind_persol2=find(xlc(end,:)>2,1,'first');

% first (unstable) solution

figure(10)
plot(mesh_refined,Per_Solutions_v(:,ind_persol1),'b',mesh_refined,xe(M+2,ind_eq39)*ones(1,length(mesh_refined)),'r');
for kk= 1:size(Per_Solutions_v,1)
    if  (Per_Solutions_v(kk,ind_persol1)<=xe(M+2,ind_eq39) && Per_Solutions_v(kk+1,ind_persol1) > xe(M+2,ind_eq39))
        ind_switch=kk;
        break;
    end
end

Profile_w_1 = [Per_Solutions_w(ind_switch+1:end,ind_persol1);Per_Solutions_w(2:ind_switch+1,ind_persol1)];
Profile_v_1 = [Per_Solutions_v(ind_switch+1:end,ind_persol1);Per_Solutions_v(2:ind_switch+1,ind_persol1)];
Profile_m_1 = [Per_Solutions_m(ind_switch+1:end,ind_persol1);Per_Solutions_m(2:ind_switch+1,ind_persol1)];

figure(3)
subplot(3,1,1)
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_w_1;Profile_w_1]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(M+1,ind_eq39);xe(M+1,ind_eq39)]);
xlabel('t','interpreter','latex');
title(['Periodic solution, w-component, p=',num2str(xlc(end,ind_persol1))...
    ', period=',num2str(xlc(end-1,ind_persol1))]);

subplot(3,1,2)
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_v_1;Profile_v_1]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(M+2,ind_eq39);xe(M+2,ind_eq39)]);
title(['Periodic solution, v-component, p=',num2str(xlc(end,ind_persol1))...
    ', period=',num2str(xlc(end-1,ind_persol1))]);


subplot(3,1,3)
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_m_1;Profile_m_1]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(M,ind_eq39);xe(M,ind_eq39)]);
xlabel('t','interpreter','latex');
title(['Periodic solution, m(t,xbar), p=',num2str(xlc(end,ind_persol1))]);

% figname=[num2str(M),'_persol1'];
% savefig(figname);

% second (stable) solution

figure(10)
plot(mesh_refined,Per_Solutions_v(:,ind_persol2),'b',mesh_refined,xe(M+2,ind_eq39)*ones(1,length(mesh_refined)),'r');
for kk= 1:size(Per_Solutions_v,1)
    if  (Per_Solutions_v(kk,ind_persol2)<=xe(M+2,ind_eq39) && Per_Solutions_v(kk+1,ind_persol2) > xe(M+2,ind_eq39))
        ind_switch=kk;
        break;
    end
end

Profile_w_2 = [Per_Solutions_w(ind_switch+1:end,ind_persol2);Per_Solutions_w(2:ind_switch+1,ind_persol2)];
Profile_v_2 = [Per_Solutions_v(ind_switch+1:end,ind_persol2);Per_Solutions_v(2:ind_switch+1,ind_persol2)];
Profile_m_2 = [Per_Solutions_m(ind_switch+1:end,ind_persol2);Per_Solutions_m(2:ind_switch+1,ind_persol2)];

figure(4)
subplot(3,1,1)
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_w_2;Profile_w_2]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(M+1,ind_eq39);xe(M+1,ind_eq39)]);
xlabel('t','interpreter','latex');
title(['Periodic solution, w-component, p=',num2str(xlc(end,ind_persol2))...
    ', period=',num2str(xlc(end-1,ind_persol2))]);

subplot(3,1,2)
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_v_2;Profile_v_2]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(M+2,ind_eq39);xe(M+2,ind_eq39)]);
title(['Periodic solution, v-component, p=',num2str(xlc(end,ind_persol2))...
    ', period=',num2str(xlc(end-1,ind_persol2))]);


subplot(3,1,3)
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_m_2;Profile_m_2]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(M,ind_eq39);xe(M,ind_eq39)]);
xlabel('t','interpreter','latex');
title(['Periodic solution, m(t,xbar), p=',num2str(xlc(end,ind_persol2))]);

% figname=[num2str(M),'_persol2'];
% savefig(figname);






function [w,x,D,q]=cheb(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% D - differentiation matrix
% q - row vector of the barycentric weights
% see Trefethen 2000, Spectral Methods in Matlab

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