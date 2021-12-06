clear all,close all,clc;
%%
% this is for SDPT3
current_dir=pwd;
cd('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
cd(current_dir);

% this is for YALMIP
addpath(genpath('C:\Users\Artun\Desktop\BACK UP 2021-june 20\Documents 2020\MATLAB toolboxes for optimization\YALMIP-master'))
yalmip('clear');
%% YALMIP dim-2 nonlinear dyn sys lyap example Linear Programming
clear all,close all,clc;yalmip('clear');
sdpvar x1 x2 
x=[x1;x2];

f1=-x2-(3/2)*x1^2-(1/2)*x1^3;
f2=3*x1-x2;

f=[f1;f2];

eps1=1e-3;

v1 = monolist([x1 x2],2,1);    sdisplay(v1)

ALPHA=sdpvar(length(v1),length(v1),'symmetric'); sdisplay(ALPHA) % 

[c,v] = coefficients(v1'*ALPHA*v1,[x]);              sdisplay([c,v])

Vx=v1'*ALPHA*v1;
[c,v] = coefficients(Vx,[x]);              sdisplay([c,v])
%% random point generation
N_random_points=1000;
neighborhood_radius=20;

[rand_vec] = generate_rnd_pts_in_ss(2,neighborhood_radius,N_random_points);

% plot(rand_vec(:,1),rand_vec(:,2),'b*'); hold on;
% plot(neighborhood_radius*cos(0:pi/100:2*pi),neighborhood_radius*sin(0:pi/100:2*pi),'r-');
% xlabel('x1');ylabel('x2');
%% V(x)>0 -> V(x)>=eps1*xi^T*xi forall i=1...N [constraints-1]
cons_vec1=[];
for ii=1:1:N_random_points
    cons1=replace(Vx,[x],[rand_vec(ii,:)]);
%     [c,v] = coefficients(cons1,[x]);              sdisplay([c,v])
    cons_vec1=[cons_vec1;coefficients(cons1,[x])>=eps1*rand_vec(ii,:)*rand_vec(ii,:)'];
end
% sdisplay(cons_vec)
size(cons_vec1)
% cons1=replace(Vx,[x],[2 2])
[c,v] = coefficients(cons1,[ALPHA]);              sdisplay([c,v])
%% Vdot computation , recall Vdor = gradV'*f(x)

Vdot=jacobian(Vx,[x])*[f];        sdisplay(Vdot)
[c,v] = coefficients(Vdot,[x]);          sdisplay([c,v])
%% V(x)dot<0 -> Vdot(x)<=-eps1*xi^T*xi forall i=1...N [constraints-1]
cons_vec2=[];
for ii=1:1:N_random_points
    cons1=replace(Vdot,[x],[rand_vec(ii,:)]);
%     [c,v] = coefficients(cons1,[x]);              sdisplay([c,v])
    cons_vec2=[cons_vec2;coefficients(cons1,[x])<=-eps1*rand_vec(ii,:)*rand_vec(ii,:)'];
end
% sdisplay(cons_vec)
size(cons_vec2)
% [c,v] = coefficients(cons1,[ALPHA]);              sdisplay([c,v])
%% YALMIP optimization problem computation

F=[];
F=[F; cons_vec1 ];
F=[F; cons_vec2 ];
F=[F; -20<=vec(ALPHA)<=20 ];

ops = sdpsettings('solver','sdpt3');
ops.dualize=[0];
ops.sos.newton=[1];
ops.sos.congruence=[1];
ops.sos.numblkdg=[1e-6];
ops.sdpt3.maxit=100;

    sol = optimize(F,[.5*norm(vec(ALPHA),1)-trace(ALPHA)],ops);
    sol.info

    ALPHA0=value(ALPHA); clean(ALPHA0,1e-5)
    min(eig(ALPHA0))
sol.info


% lambda_vec=[linspace(0.4,0.7,100)];
% nnz_vec=zeros(size(lambda_vec));
% for ii=1:1:length(lambda_vec)
%     sol = optimize(F,[lambda_vec(ii)*norm(vec(ALPHA),1)-trace(ALPHA)],ops);
%     sol.info
% 
%     ALPHA0=value(ALPHA);
%     nnz_vec(ii)=nnz(clean(ALPHA0,1e-5));
% end
% plot(lambda_vec,nnz_vec);
%%
% [min_nnz,index1]=min(nnz_vec)
%     sol = optimize(F,[lambda_vec(index1)*norm(vec(ALPHA),1)-trace(ALPHA)],ops);
%     sol.info
% 
%     ALPHA0=value(ALPHA); clean(ALPHA0,1e-5)
%% let us check if the computed V(x) satisfies the constraints
Vx0=v1'*ALPHA0*v1;
Vx0=clean(Vx0,1e-5);
[c,v] = coefficients(Vx0,[x]);              sdisplay([c,v])



Vx_values=zeros(N_random_points,1);
for ii=1:1:N_random_points
    Vx_values(ii)=replace(Vx0,[x],rand_vec(ii,:));
end
fig1=figure(1); fig1.Color=[1,1,1];
plot3(rand_vec(:,1),rand_vec(:,2),Vx_values,'b*');
xlabel('x1');ylabel('x2');zlabel('V(x)');
box on;
    
    % rand_vec(ii,:)

min(Vx_values-eps1*[rand_vec(:,1).*rand_vec(:,1)+rand_vec(:,2).*rand_vec(:,2)])
min(Vx_values)
% plot(1:1:length(Vx_values),Vx_values)
[xx1,xx2] = meshgrid(-20:0.5:20);
hold on; surf(xx1,xx2,zeros(size(xx1)));


%% let us check if the computed Vdot(x) satisfies the constraints
Vdot0=jacobian(Vx0,[x])*[f];        %sdisplay(Vdot)
Vdot0=clean(Vdot0,1e-5);
[c,v] = coefficients(Vdot0,[x]);          sdisplay([c,v])
degree(Vdot)


Vdot_values=zeros(N_random_points,1);
for ii=1:1:N_random_points
    Vdot_values(ii)=replace(Vdot0,[x],rand_vec(ii,:));
end
fig1=figure(1); fig1.Color=[1,1,1];
plot3(rand_vec(:,1),rand_vec(:,2),Vdot_values,'b*');
xlabel('x1');ylabel('x2');zlabel('Vdot');
box on;

[xx1,xx2] = meshgrid(-20:0.5:20);
hold on; surf(xx1,xx2,zeros(size(xx1)));
    % rand_vec(ii,:)

max(Vdot_values+eps1*[rand_vec(:,1).*rand_vec(:,1)+rand_vec(:,2).*rand_vec(:,2)])
max(Vdot_values)


%