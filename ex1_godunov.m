% ex1_godunov.m: 
% Minimal example of the godunov method using a background viscosity dependent 
% Riemann solver for two inviscid Burgers equations.
%
% This programm uses the matlab function riSolBurgers2.m
%
% AUTHOR:
% Valentin Pellhammer
% Department of Mathematics and Statistics,
% University of Konstanz, 78457 Konstanz
% email adress: valentin.pellhammer@uni-konstanz.de
% homepage: http://www.math.uni-konstanz.de/~pellhammer/
%
% Date: February 2020
% Last revision: February 2020

clc;close all;clear all;

% define discretisation parameters
xmin = -1;
xmax = 1;
dt = 0.001;
N = 200;
tmax = 1;
dx = (xmax-xmin)/(N-1);
x =  linspace(xmin, xmax,N);


% define model 
kappa = 0.5;
F = @(u,v) [(1/2) .* u.^2;(1/2) .* v.^2];

% Initial data %
% Riemann data (transmissive boundary cond. recomended; cf. main loop)
uL = 1;
vL = 0;
uR = -7.1;
vR = 0;
init =  [(x<0) .* uL + (x>=0) .* uR;(x<0) .* vL + (x>=0) .* vR];
% other data
% init = [cos(5.*x);sin(10.*x)];

sol=init;

%initialize plots
subplot(2,1,1)
h1 = plot(x,sol(1,:),'*');hold on;
T = title('Godunov method (press any key)');
axis([xmin xmax,-10 3]);
xlabel('x','Interpreter','latex');
ylabel('u','Interpreter','latex','Rotation',0,'FontSize',14);

subplot(2,1,2)
h2 = plot(x,sol(2,:),'*');hold on;
xlabel('$x$','Interpreter','latex','FontSize',14);
axis([xmin xmax,-3 3]);
xlabel('$x$','Interpreter','latex','FontSize',14);
ylabel('$v$','Interpreter','latex','Rotation',0,'FontSize',14);
a = annotation('textbox',[0.8,0.07,0,0],'string','$t=0$','Interpreter','latex','FontSize',14);

pause;
T.String = 'Godunov method';

for i= 0:dt:tmax
    % transmissive boundary conditions
    %solm = [sol(:,1),sol(:,1:end-1)];
    %solp = [sol(:,2:end),sol(:,end)]; 
    
    % periodic boundary conditions
    solm = [sol(:,end),sol(:,1:end-1)];
    solp = [sol(:,2:end),sol(:,1)]; 
    
    [Sp1,Sp2] = riemannSolverBurgersSquared(sol(1,:),sol(2,:),solp(1,:),solp(2,:),kappa,0);
    [Sm1,Sm2] = riemannSolverBurgersSquared(solm(1,:),solm(2,:),sol(1,:),sol(2,:),kappa,0);
    sol = sol - (dt./dx).*(F(Sp1,Sp2) - F(Sm1,Sm2));
   
    set(h1,'Ydata',sol(1,:));
    set(h2,'Ydata',sol(2,:));
    a.String = ['$t=',num2str(i),'$'
        ];

    pause(0.01);
end



