% Copyright (c) 2020 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE.txt for 
% licensing information. 
% 
% If using this code, please cite 
% Scarabel, Breda, Diekmann, Gyllenberg, Vermiglio (2020), Numerical
% bifurcation analysis of physiologically structured population models via
% pseudospectral approximation, Vietnam J Math
%
%% MC_size_Daphnia_pw
% command line instructions for Matcont continuation of the system defined
% in PS_size_Daphnia_pw

%profile on
clear; close all

M=10; 

% Piecewise:
k=2;

MM= k*M+1;

% initial value of continuation parameter
%K=0.05; mu=0.05; % chemostat
%K=0.09; mu=0.075;
K=0.1; mu=0.3; %0.1; %0.25;
mu=0.2;

aux=1;

par=[K;mu;aux;M]; % vector of initial parameters
ap1=1; % index of the bifurcation parameter in the vector 'par'
ap2=2;

% initial equilibrium point
xeq=0; yeq=K;

% tolerance
TOL=1e-3;
TestTOL=TOL;
VarTOL=TOL;

%% Continuation

handles=feval(@PS_size_Daphnia_pw); 
opt=contset; 
global cds

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
%opt=contset(opt,'MaxStepSize',0.01); % usato per M=5
opt=contset(opt,'Backward',0);

Weq = feval(handles{1},M,xeq,yeq); % initializes equilibrium vector

[x0,v0]=init_EP_EP(@PS_size_Daphnia_pw,Weq,par0,ap1);

tic
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt);
time1=toc;

xe(end,end)

% [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
% xe(end,end)

figure(2)
subplot(2,1,1)
cpl(xe,ve,se,[MM+1 MM]); hold on
xlabel('K','interpreter','latex');
title('bifurcation of S');

subplot(2,1,2)
cpl(xe,ve,se,[MM+1 MM-1]); hold on
xlabel('K','interpreter','latex');
title('bifurcation of m(t,xbar)');
%ylabel('max/min','interpreter','latex')

%% Detection of singular points
% xe,ve,se,he,fe = output of previous continuation
% par = current parameter vector
% ap1 = index of continuation parameter in vector par

% BP, branching point
for ii=1:length(se)
    if strcmp(se(ii).label,'BP')==1
        BP_index=se(ii).index;
        sBP=se(ii);
        break;
    end
end
par(ap1)=xe(end,BP_index);
BP=xe(1:end-1,BP_index);

xeBP=xe; veBP=ve; seBP=se; heBP=he; feBP=fe;
parBP=par;

%% Equilibrium continuation from BP
% xeBP = vector of variables at BP
% parBP = parameter vector at BP
% sBP = information about BP from previous continuation
disp('Starting equilibrium continuation from BP');

% set options
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Backward',0);
opt=contset(opt,'Eigenvalues',1);
opt=contset(opt,'MaxStepsize',5e-2);

    % detection of the equilibrium at K=1 or K=2
    UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P';
    opt=contset(opt,'Userfunctions',1);
    opt=contset(opt,'UserfunctionsInfo',UserInfo);
    
[x0,v0]=init_BP_EP(@PS_size_Daphnia_pw,BP,parBP,sBP,0.01);

tic
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
time2=toc

[xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)

% while ((length(se)<3) && xe(end,end)< UpperBound)
%     [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
% end

figure(2)
subplot(2,1,1)
plot(xe(end,:),xe(MM,:)); hold on
%cpl(xe,ve,se,[MM+1 MM]); hold on
xlabel('K','interpreter','latex');
title('bifurcation of S');

subplot(2,1,2)
plot(xe(end,:),xe(MM-1,:)); hold on
% cpl(xe,ve,se,[MM+1 2*M]); hold on
xlabel('K','interpreter','latex');
title('bifurcation of m(t,xbar)');
%ylabel('max/min','interpreter','latex')

% or equivalently
% plot(xe(:,end),xe(:,1))

%% Detection of singular points
% H, Hopf point

for ii=1:size(se)
    if strcmp(se(ii).label,'H ')==1
        H_index=se(ii).index;
        break;
    end
end
par(ap1)=xe(end,H_index);
H=xe(1:MM,H_index);

xeH=xe; veH=ve; seH=se; heH=he; feH=fe;
parH=par;

% P1, P2 corresponding to K=1 and K=2
for ii=1:length(se)
    if (strcmp(se(ii).label,'P')==1 && abs((xe(end,se(ii).index)-1))<0.01)
        ind_eq1=se(ii).index;
    elseif (strcmp(se(ii).label,'P')==1 && abs((xe(end,se(ii).index)-1.5))<0.01)
        ind_eq15=se(ii).index;
    elseif (strcmp(se(ii).label,'P')==1 && abs((xe(end,se(ii).index)-2))<0.01)
        ind_eq2=se(ii).index;
    end
end

% %% H continuation in two parameters
% % H = vector of variables at H
% % parH = parameter vector at H
% % ap1,ap2 = index of continuation parameters in the vector par
% display('Starting H continuation');
% 
% % set options
% TOL= 1e-6;
% opt=contset(opt,'FunTolerance',TOL); 
% opt=contset(opt,'VarTolerance',TOL);
% opt=contset(opt,'TestTolerance',TOL);
% 
%     % turn off userfunctions
%     UserInfo.name='userf1'; UserInfo.state=0;
%     opt=contset(opt,'Userfunctions',0);
%     
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'MaxNumPoints',50);
% opt=contset(opt,'Eigenvalues',1);
% opt=contset(opt,'Backward',0);
% opt=contset(opt,'MaxStepsize',1e-1);
% opt=contset(opt,'MinStepsize',1e-6);
% 
% [x0,v0]=init_H_H(@PS_size_Daphnia_pw,H,parH,[ap2 ap1]);
% [xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)
% [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
% [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
% % jj=0;
% % while (xh(MM+2,end)>0 && xh(MM+2,end)<2 && jj<5)
% %      [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
% %       jj=jj+1;
% % end
% 
% % Plot
% figure
% cpl(xh,vh,sh,[MM+1 MM+2]); hold on
% xlabel('mu'); ylabel('K');
% 
% %%
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'Eigenvalues',1);
% opt=contset(opt,'Backward',1);
% opt=contset(opt,'MaxStepsize',2e-1);
% %opt=contset(opt,'MinStepsize',1e-5);
% opt=contset(opt,'MaxNumPoints',10);
% %opt=contset(opt,'Adapt',3);
% 
% [x0,v0]=init_H_H(@PS_size_Daphnia_pw,H,parH,[ap2 ap1]);
% [xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(MM+1,end)
% [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
% jj=0;
% while (xhb(MM+1,end)>0 && xhb(MM+1,end)<0.32 && xhb(MM+2,end)>0 && jj<5)
%      [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
%       jj=jj+1;
% end
% cpl(xhb,vhb,shb,[MM+1 MM+2]); hold on
% %title(['Hopf curve, M=',num2str(M)]);
% 
% %figname=[num2str(M),'_hopf'];
% %savefig(figname);
% 
% %%
% ap_aux=3;
% [xb0,vb0]=init_BP_BP(@PS_size_Daphnia_pw,xeBP(1:end-1,BP_index),parBP,[ap2 ap1 ap_aux],ap2);
% opt=contset(opt,'MaxNumPoints',100);
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'Backward',0);
% [xb1,vb1,sb1,hb1,fb1]=cont(@branchpoint,xb0,vb0,opt); xb1(end-2,end)
% [xb1,vb1,sb1,hb1,fb1]=cont(xb1,vb1,sb1,hb1,fb1,cds); xb1(end-2,end)
% jj=1;
% while (xb1(MM+1,end)>=0 && xb1(MM+1,end)<=0.345 && xb1(MM+2,end)>=0 && jj<10)
%     [xb1,vb1,sb1,hb1,fb1]=cont(xb1,vb1,sb1,hb1,fb1,cds); xb1(MM+1,end)
%     jj=jj+1;
% end
% %figure
% cpl(xb1,vb1,sb1,[MM+1 MM+2]); hold on;
% 
% %%
% ap_aux=3;
% opt=contset(opt,'Backward',1);
% [xb0,vb0]=init_BP_BP(@PS_size_Daphnia_pw,xeBP(1:end-1,BP_index),parBP,[ap2 ap1 ap_aux],ap2);
% [xb2,vb2,sb2,hb2,fb2]=cont(@branchpoint,xb0,vb0,opt);
% jj=1;
% while (xb2(MM+1,end)>=0 && xb2(MM+1,end)<=0.34 && xb2(MM+2,end)>=0 && jj<5)
%      [xb2,vb2,sb2,hb2,fb2]=cont(xb2,vb2,sb2,hb2,fb2,cds);
%      jj=jj+1;
% end
% cpl(xb2,vb2,sb2,[MM+1 MM+2]);
% xlabel('mu'); ylabel('K');
% %title(['transcritical curve, M=',num2str(M)]);
% %figname=[num2str(M),'_transcritical'];
% 
% title(['Stability regions, piecewise, M=',num2str(M)]);
% axis([0 0.5 0 12])
% figname=[num2str(M),'_regions_piecewise'];
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
% opt=contset(opt,'MinStepsize',1e-6);

TOL=1e-6;
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);

opt=contset(opt,'MaxNumPoints',50);
opt=contset(opt,'MaxStepsize',1);
% opt=contset(opt,'MinStepsize',1e-16);
%opt=contset(opt,'InitStepsize',1e-2);
%opt=contset(opt,'Adapt',0);
%opt=contset(opt,'MaxNumPoints',50);

%opt=contset(opt,'Increment',1e-6); %default 1e-5
%opt=contset(opt,'MaxCorrIters',100); %default 1e-5
%opt=contset(opt,'MaxNewtonIters',10); %default 1e-5
%opt=contset(opt,'Increment',1e-10); %default 1e-5

    % detection of the equilibrium at K=1 or K=2
    UserInfo.name='userf'; UserInfo.state=1; UserInfo.label='P';
    opt=contset(opt,'Userfunctions',1);
    opt=contset(opt,'UserfunctionsInfo',UserInfo);

ntst=40; % number of intervals
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_size_Daphnia_pw,H,parH,ap1,1e-3,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
[xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)

jj=0;
while (max(xlc(end,:))<2 && jj<5)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
    jj=jj+1;
end

%% Plot max and min periodic solutions
    % re-interpolation for a smoother plot
    
    mesh_refined=linspace(0,1,100);
    Per_Solutions_m = zeros(length(mesh_refined),size(xlc,2));
    Per_Solutions_S = zeros(length(mesh_refined),size(xlc,2));
    
    xb=0.8; xA=2.5; xm=6; 
    [~,Nodes,~,BaryWeights]=cheb(M,xb,xm);

  xlc_j = zeros(ntst+1,1);
  Per_Solutions_j = zeros(length(mesh_refined),size(xlc,2)); % juveniles
  Per_Solutions_a = zeros(length(mesh_refined),size(xlc,2)); % adults

    for ind_persol=1:size(xlc,2)

        Per_Solutions_m(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(MM-1:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        Per_Solutions_S(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(MM:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        Per_Solutions_j(:,ind_persol) = interp1(flc(1:ntst+1,ind_persol),xlc(M:MM*ncol:((ntst*ncol+1)*MM),ind_persol),mesh_refined,'spline');
        Per_Solutions_a(:,ind_persol) = Per_Solutions_m(:,ind_persol)-Per_Solutions_j(:,ind_persol);
    
    end
    
upperbound_m=max(Per_Solutions_m);
lowerbound_m=min(Per_Solutions_m);
upperbound_S=max(Per_Solutions_S);
lowerbound_S=min(Per_Solutions_S);

index=size(xlc,2);

figure(2)
subplot(2,1,1)
plot(xe(end,:),xe(MM,:),'b'); hold on
for jj=1:size(se,1)
    plot(xe(end,se(jj).index),xe(MM,se(jj).index),'.r');
end
%cpl(xe,ve,se,[MM+1 MM]); hold on
plot(xlc(end,1:index),upperbound_S(1:index),'g',xlc(end,1:index),lowerbound_S(1:index),'g');
xlabel('K','interpreter','latex');
title('bifurcation of S');
axis([0 2 0 3])

subplot(2,1,2)
plot(xe(end,:),xe(MM-1,:),'b'); hold on
for jj=1:size(se,1)
    plot(xe(end,se(jj).index),xe(MM-1,se(jj).index),'.r');
end
%cpl(xe,ve,se,[MM+1 M]); hold on
plot(xlc(end,1:index),upperbound_m(1:index),'g',xlc(end,1:index),lowerbound_m(1:index),'g');
xlabel('K','interpreter','latex');
title('bifurcation of m(t,xbar)');
axis([0 2 0 0.05])
%ylabel('max/min','interpreter','latex')



% figname=[num2str(M),'_bif_Daphnia'];
% savefig(figname);

% %% Plot periodic profiles in normalized period
% 
% % P1, P2 corresponding to K=1 and K=2
% for ii=1:length(slc)
%     if strcmp(slc(ii).label,'P1')==1
%         ind_persol1=slc(ii).index;
%     elseif strcmp(slc(ii).label,'P2')==1
%         ind_persol2=slc(ii).index;
%     end
% end
% 
%
% % first solution
% ind_persol1=find(xlc(end,:)>1,1,'first');
% figure(3)
% subplot(2,1,1)
% plot([mesh_refined,1+mesh_refined],[Per_Solutions_m(:,ind_persol1);Per_Solutions_m(:,ind_persol1)]);
% title(['Periodic solution, m(t,xbar), K=',num2str(xlc(end,ind_persol1)),...
%     ', period=',num2str(xlc(end-1,ind_persol1))]);
% 
% subplot(2,1,2)
% plot([mesh_refined,1+mesh_refined],[Per_Solutions_S(:,ind_persol1);Per_Solutions_S(:,ind_persol1)]);
% xlabel('K','interpreter','latex');
% title(['Periodic solution, S-component, K=',num2str(xlc(end,ind_persol1))]);
% figname=[num2str(M),'_persol1'];
% savefig(figname);
% 
% % second solution
% ind_persol2=find(xlc(end,:)>2,1,'first');
% figure(4)
% subplot(2,1,1)
% plot([mesh_refined,1+mesh_refined],[Per_Solutions_m(:,ind_persol2);Per_Solutions_m(:,ind_persol2)]);
% title(['Periodic solution, m(t,xbar), K=',num2str(xlc(end,ind_persol2)),...
%     ', period=',num2str(xlc(end-1,ind_persol2))]);
% 
% subplot(2,1,2)
% plot([mesh_refined,1+mesh_refined],[Per_Solutions_S(:,ind_persol2);Per_Solutions_S(:,ind_persol2)]);
% xlabel('K','interpreter','latex');
% title(['Periodic solution, S-component, K=',num2str(xlc(end,ind_persol2))]);
% figname=[num2str(M),'_persol2'];
% savefig(figname);

% figure
% plot(flc(ntst+2:end,ind_persol1),'.b')
% hold on
% plot(flc(ntst+2:end,ind_persol2),'.r')
% title('multipliers')
% legend(['K=',num2str(xlc(end,ind_persol1))],['K=',num2str(xlc(end,ind_persol2))])


%% Plot periodic profiles in real period

% P1, P2 corresponding to K=1 and K=2
for ii=1:length(slc)
    if (strcmp(slc(ii).label,'P')==1 && abs(xlc(end,slc(ii).index)-1.5)<0.05)
        ind_persol1=slc(ii).index;
    elseif (strcmp(slc(ii).label,'P')==1 && abs(xlc(end,slc(ii).index)-2)<0.05)
        ind_persol2=slc(ii).index;
    end
end
% 
% ind_persol1=find(xlc(end,:)>1.5,1,'first');
% ind_persol2=find(xlc(end,:)>2,1,'first');

% first solution

figure(10)
plot(mesh_refined,Per_Solutions_S(:,ind_persol1),'b',mesh_refined,xe(MM,ind_eq15)*ones(1,length(mesh_refined)),'r');
for kk= 1:size(Per_Solutions_S,1)
    if  (Per_Solutions_S(kk,ind_persol1)<=xe(MM,ind_eq15) && Per_Solutions_S(kk+1,ind_persol1) > xe(MM,ind_eq15))
        ind_switch=kk;
        break;
    end
end

Profile_S_1 = [Per_Solutions_S(ind_switch+1:end,ind_persol1);Per_Solutions_S(2:ind_switch+1,ind_persol1)];
Profile_m_1 = [Per_Solutions_m(ind_switch+1:end,ind_persol1);Per_Solutions_m(2:ind_switch+1,ind_persol1)];
Profile_j_1 = [Per_Solutions_j(ind_switch+1:end,ind_persol1);Per_Solutions_j(2:ind_switch+1,ind_persol1)];
Profile_a_1 = [Per_Solutions_a(ind_switch+1:end,ind_persol1);Per_Solutions_a(2:ind_switch+1,ind_persol1)];

figure(3)
subplot(3,1,1)
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_S_1;Profile_S_1]); hold on
plot([0;2*xlc(end-1,ind_persol1)],[xe(MM,ind_eq15);xe(MM,ind_eq15)]);
xlabel('t','interpreter','latex');
title(['Periodic solution, S-component, K=',num2str(xlc(end,ind_persol1))...
    ', period=',num2str(xlc(end-1,ind_persol1))]);

subplot(3,1,2)
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_m_1;Profile_m_1]); hold on
plot([0;2*xlc(end-1,ind_persol1)],[xe(MM-1,ind_eq15);xe(MM-1,ind_eq15)]);
title(['Periodic solution, m(t,xbar), K=',num2str(xlc(end,ind_persol1))]);

subplot(3,1,3)
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_j_1;Profile_j_1]); hold on
plot(xlc(end-1,ind_persol1)*[mesh_refined,1+mesh_refined],[Profile_a_1;Profile_a_1]);
xlabel('t','interpreter','latex');
legend('juv','adu');
title(['Periodic solution, juveniles and adults, K=',num2str(xlc(end,ind_persol1))]);

% figname=[num2str(M),'_persol1'];
% savefig(figname);

% second solution

for kk=1:size(Per_Solutions_S,1)-1
    if  (Per_Solutions_S(kk,ind_persol2)<xe(MM,ind_eq2)) & (Per_Solutions_S(kk+1,ind_persol2)>xe(MM,ind_eq2))
        ind_switch2=kk;
        break;
    end
end

Profile_S_2 = [Per_Solutions_S(ind_switch2+1:end,ind_persol2);Per_Solutions_S(2:ind_switch2+1,ind_persol2)];
Profile_m_2 = [Per_Solutions_m(ind_switch2+1:end,ind_persol2);Per_Solutions_m(2:ind_switch2+1,ind_persol2)];
Profile_j_2 = [Per_Solutions_j(ind_switch2+1:end,ind_persol2);Per_Solutions_j(2:ind_switch2+1,ind_persol2)];
Profile_a_2 = [Per_Solutions_a(ind_switch2+1:end,ind_persol2);Per_Solutions_a(2:ind_switch2+1,ind_persol2)];


figure(4)
subplot(3,1,1)
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_S_2;Profile_S_2]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(MM,ind_eq2);xe(MM,ind_eq2)]);
xlabel('K','interpreter','latex');
title(['Periodic solution, S-component, K=',num2str(xlc(end,ind_persol2)),...
', period=',num2str(xlc(end-1,ind_persol2))]);

subplot(3,1,2)
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_m_2;Profile_m_2]); hold on
plot([0;2*xlc(end-1,ind_persol2)],[xe(MM-1,ind_eq2);xe(MM-1,ind_eq2)]);
title(['Periodic solution, m(t,xbar), K=',num2str(xlc(end,ind_persol2))]);

subplot(3,1,3)
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_j_2;Profile_j_2]); hold on
plot(xlc(end-1,ind_persol2)*[mesh_refined,1+mesh_refined],[Profile_a_2;Profile_a_2]);
legend('juv','adu');
title(['Periodic solution, juveniles and adults, K=',num2str(xlc(end,ind_persol2))]);

% figname=[num2str(M),'_persol2'];
% savefig(figname);



% filename=[num2str(M),'_bif_Daphnia'];
% save(filename);


% %% Plot stationary distributions
% 
% % ind_eq1=find(xe(end,:)>1,1,'first');
% % ind_eq2=find(xe(end,:)>2,1,'first');
% 
%     xb=0.8; xA=2.5; xm=6; 
%     [QuadWeights,Nodes,DD,BaryWeights]=cheb(M,xb,xm);
% 
%     dMDM(1:M,1:M+1)= DD(2:end,:);
% 	DM = DD(2:end,2:end);
%     
%     mesh_eq=Nodes; %linspace(xb,xm,21);
% 
% % prof_n1=interpoly(mesh_eq,Nodes,DD(1:end,2:end)*xe(1:M,ind_eq1),BaryWeights);
% % prof_n2=interpoly(mesh_eq,Nodes,DD(1:end,2:end)*xe(1:M,ind_eq2),BaryWeights);
% % figure
% % plot(mesh_eq,prof_n1,mesh_eq,prof_n2);
% % title('Stationary distribution n^*(x)');
% % legend(['K=',num2str(xe(end,ind_eq1))],['K=',num2str(xe(end,ind_eq2))])
% 
% prof_m1=interpoly(mesh_eq,Nodes,[0;xe(1:M,ind_eq1)],BaryWeights);
% prof_m2=interpoly(mesh_eq,Nodes,[0;xe(1:M,ind_eq2)],BaryWeights);
% 
% figure
% subplot(1,2,1)
% plot(mesh_eq,prof_m1);
% xlabel('x');
% title(['Stationary distribution m(x), K=', num2str(xe(end,ind_eq1))]);
% %legend(['K=',num2str(xe(end,ind_eq1))],['K=',num2str(xe(end,ind_eq2))])
% 
% subplot(1,2,2)
% plot(mesh_eq,prof_m2);
% xlabel('x');
% title(['Stationary distribution m(x), K=', num2str(xe(end,ind_eq2))]);
% 
% figname=[num2str(M),'_eq_distr'];
% savefig(figname);
% 
% %% Analytic profiles:
% mesh_eq=linspace(xb,xm,101);
% [m_bar_1,s_bar_1] = Equil_Daphnia(1,mu,mesh_eq);
% [m_bar_2,s_bar_2] = Equil_Daphnia(2,mu,mesh_eq);
% 
% figure
% subplot(1,2,1)
% plot(mesh_eq,m_bar_1)
% title('K=1')
% subplot(1,2,2)
% plot(mesh_eq,m_bar_2)
% title('K=2')
% figname=['0_eq_distr_analytic'];
% savefig(figname);

function [w,x,D,q]=cheb(N,a,b)
% Output:
% x - N+1 Chebyshev nodes on [a,b] (x_0=a, x_N=b),
% w - weights of the quadrature formula in [a,b],
% D - differentiation matrix
% q - row vector of the barycentric weights
% see Trefethen 2000, Spectral methods in Matlab

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