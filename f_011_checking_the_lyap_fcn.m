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
%% 
clear all,close all,clc;yalmip('clear');

tspan = [0:0.1:10];
x_init = [1;1];

odefcn1=@(t,x,params) [-x(2)-(3/2)*x(1)^2-(1/2)*x(1)^3;3*x(1)-x(2)];
params=1;
[t,x] = ode45(@(t,x) odefcn1(t,x,params), tspan, x_init);

% x=[x1;x2]; 
% f1=-x2-(3/2)*x1^2-(1/2)*x1^3;
% f2=3*x1-x2;
% f=[f1;f2];
fig1=figure(1);fig1.Color=[1,1,1];
% plot3(t,x(:,1),x(:,2),'linewidth',[2]);
% xlabel('t');ylabel('x_1(t)');zlabel('x_2(t)');
% % plot(x(:,1),x(:,2),'linewidth',[2]);
% % xlabel('x_1(t)');ylabel('x_2(t)');
% % box on; grid on;
%% let us generate 100 random initial-points in state space [inside the Ball w/ radius=10]
% plot the "ball" [in 2-dim, the ball is a circle]
rad_ball=10;
fig1=figure(1);fig1.Color=[1,1,1];
plot(rad_ball*cos(-pi:0.1:pi)',rad_ball*sin(-pi:0.1:pi)','r-','linewidth',[2]);hold on;
% now generate random points in the ball and simulate them as
% init-conditions

Npoints=100;
% [XX] = generate_rnd_pts_in_ss(Ndim,r1,Npoints)
[rand_vec] = generate_rnd_pts_in_ss([2],rad_ball,Npoints);

for ii=1:1:Npoints
    tempx_init=rand_vec(ii,:);
    tspan = [0:0.05:5];
    [t,x] = ode45(@(t,x) odefcn1(t,x,params), tspan, tempx_init);
    plot(x(:,1),x(:,2),'linewidth',[2]); hold on;
end
xlabel('x_1(t)');ylabel('x_2(t)');
box on; grid on;
%% PLOT V(x) on x1-x2 plane
% f = @(x,y) -(x-10).^2 -(y-10).^2;
% [X,Y] = meshgrid(-10:10);
% figure(1)
% meshc(X, Y, f(X,Y))
% grid on

[XX1,XX2] = meshgrid(-10:.1:10);
% meshc(XX1, XX2, fcn_Vx(XX1,XX2))
contour(XX1,XX2,fcn_Vx(XX1,XX2),500)


% our V(x)
%     {'19.9999999958' }    {'x2^2'     }
%     {'1.84109637307' }    {'x2^4'     }
%     {'-5.18182970736'}    {'x1*x2^3'  }
%     {'19.9999999966' }    {'x1^2'     }
%     {'44.4794213698' }    {'x1^2*x2'  }
%     {'19.9999999699' }    {'x1^2*x2^2'}
%     {'-2.54758079843'}    {'x1^3*x2'  }
%     {'19.9999999932' }    {'x1^4'     }
%
function [Vx_value]=fcn_Vx(x1,x2)
%     x1=x(1);    x2=x(2);
    term1=19.9999999958*x2.^2;
    term2=1.84109637307*x2.^4;

    term3=-5.18182970736*x1.*x2.^3;
    term4=19.9999999966*x1.^2;

    term5=44.4794213698*x1.^2.*x2;
    term6=19.9999999699*x1.^2.*x2.^2;

    term7=-2.54758079843*x1.^3.*x2;
    term8=19.9999999932*x1.^4;
    Vx_value=term1+term2+term3+term4+term5...
        +term6+term7+term8;
end



% function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
% end












%






%