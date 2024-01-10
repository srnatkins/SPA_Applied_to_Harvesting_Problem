%% File Description
% This file generates figures relating to our results when using SPA to
% improve the switches to h_{4,sing} that serves as an alternative
% harvesting strategy to the following harvesting problem. 

%Problem accessed from the following paper
%Title: Marine reserves and optimal harvesting
%author: Michael G. Neubert 
%year: 2003 
%journal: Ecology letters
%PROBLEM
%max_h  \int^L_0 h(x)u(x) dx 
%               0\leq h(x)\leq h_max
%               u''(x) = -u(1-u)+h(x)*u
%               u(0) =u(L)=0

% We use y1(x)= u(x) and y2(x) = u'(x). 

% The script file already inputs the optimal switches s^*, the control form
% h_{4,sing}, and the values for y2(0) and lambda1(L). All of which were 
% numerically obtained in another script file- SPA4.m 
% 
% The file outputs two figures
% Figure 1 displays h_{4,sing}^* and the resulting state solution
% y1(x)=u(x)
% Figure 2 displays the resulting switching function and vertical lines
% indicating the locations of the switches. 

% The SuiteOPT software is not needed for this file. 
%-------------------------------------------------------------------------
function main
clf 
clear all
%% ------ Initialize parameters and variables -------%%
global L hmax h0 h1 h2 h3 h4 h N y20 lam1L ds dxstep errtol dxlam
L = 10;                            % length of habitat
hmax = 1;                          % maximum harvesting
h0 = hmax;                         % h(x) = hmax on (s0,s1) where s0 = 0      
h1 = 0;                            % h(x) = 0    on (s1,s2)
h2 = 0.5;                          % h(x) = 1/2  on (s2,s3) <- singular arc
h3 = 0;                            % h(x) = 0    on (s3,s4) 
h4 = hmax;                         % h(x) = hmax on (s4,s5) where s5 = L

h = [h0; h1; h2; h3; h4];      
N = length(h);                          
y20 =0.191999750652189;                 % y2(0)                 
y20save = y20;
lam1L =  1.302113723433043;             % lambda1(L)     
ds=L/500;
dxstep =1.e-3;
errtol =1.e-10;
dxlam = dxstep;
% The switches obtained from using SPA
s1 = 1.751269679932279;                 % sopt(1)             
s2 = 2.845228609459884;                 % sopt(2)
s3 = 7.154771390546768;                 % sopt(3)
s4 = 8.248730320073570;                 % sopt(4)
s = [s1;s2;s3;s4];
xN = 200;
[y20up,X,Y,yS] = y20update(s);
[X,Y,yS,y1L] = ysolve2(y20up,s,dxstep);
[X,Y,yS,y1L] = ysolve2(y20up,s,dxstep/4);
S =[0;s;L];
for k = 1:N
    xtest{k} = X{k};
    Ytest = Y{k};
    y1test{k} = Ytest(:,1);
    htest{k} = h(k)*ones(length(Ytest(:,1)),1);
end

for k =1:N
    xstep{k} = linspace(S(k),S(k+1),xN);
    Ynew = Y{k};
    y1{k} = Ynew(:,1);
    uint{k} = spline(X{k},y1{k},xstep{k});
end

y1value = gety1(s);
[lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(s);
[Xlam, Lam, lam20, lamS] = lamsol2(y1value, s, dxlam);
[Xlam, Lam, lam20, lamS] = lamsol2(y1value, s, dxlam/4);
for k = 1:N
    Lamnew = Lam{k};
    lam2{k} = Lamnew(:,2);
    lam2int{k} = spline(Xlam{k},lam2{k},xstep{k});
end
k=1;
i=1;
for k = 1: N
    lambda2 = lam2int{k};
    u = uint{k};
    for i=1:xN
        switching(i) = (lambda2(i)-1)*u(i);
    end
    psi{k} = switching;
end
% Plotting the Switching function
figure(1) 
plot(xstep{1},psi{1}, 'k', 'Linewidth',2);
hold on 
plot(xstep{2},psi{2},'k', 'Linewidth',2);
hold on 
plot(xstep{3},psi{3},'k', 'Linewidth',2);
hold on 
plot(xstep{4},psi{4},'k', 'Linewidth',2);
hold on 
plot(xstep{5}, psi{5},'k', 'Linewidth',2);
hold on 
yline(0,'g--');
hold on 
xline(s(1),'--r');
hold on 
xline(s(2),'--r');
hold on 
xline(s(3),'--r');
hold on 
xline(s(4),'--r');
xlabel('Habitat Location $x$','Interpreter','latex','FontSize', 25)
ylabel('$\psi(x)$','Interpreter', 'latex', 'FontSize',25);
title('Switching Function ($h_{4,sing}^*$)', 'Interpreter', 'latex', 'FontSize', 25);
ax = gca;
ax.FontSize = 20;
savefig('SPA4PsiPlot.fig');

%Control and State
figure(2) 
plot(xtest{1},y1test{1},'b','Linewidth',2)
hold on
plot(xtest{1},htest{1},'r','Linewidth',2)
hold on
plot(xtest{2},y1test{2}, 'b','Linewidth',2)
hold on
plot(xtest{2},htest{2}, 'r','Linewidth',2)
hold on
plot(xtest{3},y1test{3},'b','Linewidth',2)
hold on 
plot(xtest{3},htest{3},'r','Linewidth',2)
hold on 
plot(xtest{4},y1test{4},'b','Linewidth',2)
hold on 
plot(xtest{4},htest{4},'r','Linewidth',2)
hold on 
plot(xtest{5}, y1test{5},'b','Linewidth',2)
hold on 
plot(xtest{5}, htest{5},'r','Linewidth',2)
hold on
% The last bit is plotting the jumps of the control
plot([s(1),s(1)],[0,1],'r', 'Linewidth', 2);
hold on 
plot([s(2),s(2)],[0,0.5],'r', 'Linewidth', 2);
hold on 
plot([s(3),s(3)],[0,0.5],'r', 'Linewidth', 2);
hold on 
plot([s(4),s(4)],[0,1],'r', 'Linewidth', 2);
ylim([-.001,1.001])
xlabel('Habitat Location $x$','Interpreter','latex','FontSize', 25)
title('$h_{4,sing}^*$ and $y_1^*$', 'Interpreter', 'latex', 'FontSize', 25);
ax = gca;
ax.FontSize = 20;
%savefig('SPA4Plot.fig');





%% Functions needed for file
%% Functions for updating the initial value y20
%---------------dify - setting up ode for state equations ------------%
function  dy = dify(x,y,h)
    dy(1) = y(2);
    dy(2) = -y(1)*(1-y(1))+h*y(1);
    dy(3) = y(4);
    dy(4) = (2*y(1)-1+h)*y(3);
    dy(5) = -h*y(1);          %RHS is integrand of cost functional 
    dy = dy';
end
%------------------Using ode45 for solving dify--------------%
% we will use a for loop and cell array commands to save ode45 solution of
% dify over the intervals [S(k),S(k+1)] where S = {0,s1,...,s_N-1,L}
function [X,Y,yS, y1L, y3L] = ysolve(Y20, s)
    S = [0; s; L];
    k = 1;                      %loop counter
    yS = zeros(length(S),5);    %initializing yS- evaluates y at points in S
    yS(1,:) = [0,Y20,0,1,0];      %[y1(0),y2(0),y3(0),y4(0)] 
    intsize = zeros(N,1);       % numb. of mesh pts for solving ode over each interval
    for k = 1 : N
        options = odeset('RelTol',1.e-10,'AbsTol',errtol,'MaxStep', 1.e-2);
        ODE = @(x,y)dify(x,y,h(k));
        [X{k},Y{k}] = ode45(ODE,[S(k),S(k+1)], yS(k,:), options);
        %xnew = X{k};
        %intsize(k) = length(xnew);
        ynew = Y{k};
        yS(k+1,:) =ynew(end,:); %[y1(sk),y2(sk),y3(sk),y4(sk)];
    end
    y1L = yS(length(S),1);
    y3L = yS(length(S),3);
end
%==================NEWTON's METHOD FOR UPDATING y20=====================%
function[y20up,X,Y,yS] = y20update(s)
    tol = 1.e-12;
    maxits = 500;
    k = 1;
    oldy20 = y20;
    %fprintf (1,'newton guess: %e s: %e %e %e %e\n',y20, s) ;
    while (k<= maxits)
        [X,Y, oldyS, oldy1L, oldy3L] = ysolve(oldy20,s);
        newton = -oldy1L/oldy3L;
        if ( abs(newton) >= 2.e-4 )
            newton = 2.e-4*newton/abs(newton) ;
        end
        newy20 = oldy20+newton ;
%fprintf(1,'num: %e den: %e newton: %e y20: %e\n',oldy1L, oldy3L, newton, newy20);
        if abs(newton)< tol
            y20up = newy20;
            [X,Y,yS, y1L, y3L] = ysolve(y20up,s);
            y20 = y20up;
            break
        end
        k=k+1;
        oldy20 = newy20;
    end
    if k == maxits
        y20up = 'Maximum number of iterations exceeded';
        disp(y20up);
    end
end
%% Function for Spline interpolation 
%------------ ode system of just  u and v 
function dy2 = dify2(x,y,h)
    dy2(1) = y(2);
    dy2(2) = -y(1)*(1-y(1))+h*y(1);
    dy2 = dy2';
end
function [X,Y,yS,y1L] =ysolve2(Y20,s,dx)
    S = [0;s;L];
    k = 1;
    yS = zeros(length(S),2);
    yS(1,:) = [0,Y20];
    for k = 1:N
        options = odeset('RelTol',1.e-10,'AbsTol',errtol,'MaxStep',dx);
        ODE = @(x,y)dify2(x,y,h(k));
        [X{k},Y{k}] =ode45(ODE,[S(k),S(k+1)],yS(k,:),options);
        ynew = Y{k};
        yS(k+1,:)=ynew(end,:);
    end
    y1L = yS(length(S),1);
end 
%-----------Interpolating y1 to use for solving costate eqns-----------%
% function will use spline commands to interpolate
function y1int = y1interpolate(X,Y)
    k = 1;
    for k = 1:N
        Ynew = Y{k};
        y1{k} = Ynew(:,1);
        y1int{k} = spline(X{k},y1{k});
    end
end
% ---------gety1 interpolating function
% this function will call all of the preceeding functions to get a y1
% interpolation
function y1int = gety1(s)
    [y20up,X,Y,yS] = y20update(s);
    [X,Y,yS,y1L] = ysolve2(y20up,s,dxstep);
    [X,Y,yS,y1L] = ysolve2(y20up,s,dxstep/4);
    y1int = y1interpolate(X,Y);
end
%% Functions for the Costate
%-------costate - function sets up ode system for costate------- %%
function dlambda = costate(x,lambda, y1, h)
    dlambda(1) = h - lambda(2)*(-1+2*ppval(y1,x)+h);
    dlambda(2) = -lambda(1);
    dlambda(3) = (1-2*ppval(y1,x)-h)*lambda(4);
    dlambda(4) = -lambda(3); 
    dlambda = dlambda'; 
end
% Solving for costate equations
function [Xlam, Lam, lam20, lam40, lamS] = lamsol(Lam1L,s)
    y1int = gety1(s);
    S= [0;s;L];
    lamS = zeros(length(S),4);
    lamS(length(S),:)= [Lam1L,0,1,0];
    k = N+1;
    for k = N+1: -1: 2
        options = odeset('RelTol',1.e-10,'AbsTol',errtol,'MaxStep',1.e-2);
        ODE = @(x,lambda)costate(x,lambda,y1int{k-1},h(k-1));
        [Xlam{k-1},Lam{k-1}]=ode45(ODE, [S(k),S(k-1)], lamS(k,:),options);
        newlam = Lam{k-1};
        lamS(k-1,:) = newlam(end,:);
    end
    lam20 = lamS(1,2);
    lam40 = lamS(1,4);
end
%------Using Newton's Method for Updating lam1L-------%
function [lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(s)
    maxits = 5;
    tol = 1.e-12;
    j = 1;
    oldlam1L = lam1L;
    %Y20 = y20;
    [y20up, X, Y,yS] = y20update(s);
    while (j<maxits)
        [Xlam, Lam, oldlam20, oldlam40, lamS] = lamsol(oldlam1L,s);
        newton = -oldlam20/oldlam40;
        newlam1L = oldlam1L + newton;
        if abs(newton)<tol
            lam1Lup = newlam1L;
            [Xlam, Lam, lam20, lam40, lamS] = lamsol(lam1Lup,s);
            lam1L = lam1Lup;
            break
        end
        j = j+1;
        oldlam1L = newlam1L;
    end
    if j == maxits
        lam1Lup = 'maximum number of iterations exceeded';
        disp(lam1Lup)
    end
end
%%Functions for Computing Cost, Hamiltonian and gradient of cost
%------------- Computing Cost
function cost = objective(s)
    [y20up, X, Y,yS] = y20update(s);
    cost = yS(end,end);
end
%---------Computing Hamiltonian---------%
function H = Ham(y1, y2, lam1, lam2,h)
    lam3 = 1;
    H = lam1*y2+lam2*(-y1*(1-y1))+lam2*h*y1+lam3*(-h*y1);
end
%----------Computing Gradient of Cost wrt switching points----------%
function gradient = grad(s)
    [y20up, X, Y,yS] = y20update(s);
    [lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(s);
    lam1s = lamS(:,1);
    lam2s = lamS(:,2);
    y1s = yS(:,1);
    y2s = yS(:,2);
    i=1;
    for i = 1:length(s)
         gradient(i) = Ham(y1s(i+1),y2s(i+1),lam1s(i+1),lam2s(i+1),h(i))-Ham(y1s(i+1),y2s(i+1),lam1s(i+1),lam2s(i+1),h(i+1));
    end
end
%% Function for interpolating lam2
%--------ode system
function dlambda2 = costate2(x, lambda, y1,h)
    dlambda2(1) = h - lambda(2)*(-1+2*ppval(y1,x)+h);
    dlambda2(2) = -lambda(1);
    dlambda2 = dlambda2';
end
% solving ode system 
function [Xlam, Lam, lam20, lamS] = lamsol2(y1int, s, deltax)
    %Should have updated y20 and lambda1L and interpolate y1int
    dx = deltax;
    S =[0;s;L];
    lamS = zeros(length(S),2);
    lamS(length(S),:) = [lam1L,0];
    k = N+1;
    for k = N+1: -1:2
        options = odeset('RelTol',1.e-10,'AbsTol',errtol,'MaxStep',dx);
        ODE = @(x,lambda) costate2(x, lambda, y1int{k-1}, h(k-1));
        [Xlam{k-1},Lam{k-1}]= ode45(ODE, [S(k),S(k-1)],lamS(k,:),options);
        newlam = Lam{k-1};
        lamS(k-1,:) = newlam(end,:);
    end
    lam20 = lamS(1,2);
end
% % interpolating y1
% function lam2int = lam2interpolate(Xlam,Lam)
%     k=1; 
%     for k = 1:N
%         Lamnew = Lam{k};
%         lam2{k} = Lamnew(:,2);
%         lam2int{k} = spline(Xlam{k},lam2{k});
%     end
% end
% now we need to getlam2 by calling all of these functions
% function lam2int = getlam2(s)
%     y1int = gety1(s);
%     [lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(s);
%     [Xlam, Lam, lam20, lamS] = lamsol2(y1int, s, dxlam);
%     [Xlam, Lam, lam20, lamS] = lamsol2(y1int, s, dxlam/2);
%     lam2int = lam2interpolate(Xlam,Lam);
% end
% %% Computing the Switching function
% function [x,psi] = switching(s,y1int,lam2int)
%     S = [0;s;L];
%     i = 1;
%     k =1;
%     xN = 250;
%     for k =1:N
%         X{k} = linspace(S(k),S(k+1),xN); %partitioning [S(k),S(k+1)]
%         xuse = X{k};
%         for i =1:xN
%             PSIuse(i) = ppval(y1int{k},xuse(i));
%         end
%         PSI{k} = PSIuse;
%     end
%     x = vertcat(X{:});
%     psi = vertcat(PSI{:});
% end

end