%%---------------------------Intro.------------------------- 
% This script file uses PASA to solve for Neubert's harvesting problem via 
%a switch point algorithm in the case where the habitat length L=10 and
% the parameter hmax=1. In this file we assume the control form has 
% 4 switches.

% The initial guess of the switches, y2(0), and lambda1(L) came from using
% total variation regularization to solve the harvesting problem when the
% regularization tuning parameter is set to 10^{-3}. The values are
% computed in harvest3.m.

% You will need to download SuiteOPT software in order to run this file
% This software can be accessed on
% https://people.clas.ufl.edu/hager/software/ 
% In our experiments we used the SuiteOPT version 2.0.3. You can use the 
% same link but click software archive to access this older version of the
% software. 

% OUTPUTS
% PASA stats for each iteration of solving the SPA subproblem
%   -We also print each update of the switches and the box constraints for
%    the switches.
% sopt          - vector of the optimal switches 
% y20           - updated value of y2(0)=u'(0)
% lambda1L      - update value of lambda1(L)=lambda_u(L)
% cost          - value of the cost functional

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


% Note in our file y1(x)=u(x), y2(x)=u'(x). Additionally for our state
% equations we introduce two other differential equations which are needed
% when using Newton's method. 
%y1'(x) = y2
%y2'(x) = -y1(1-y1)+h(x)y1
%y3'(x)=y4
%y4'(x) = [2*y1-1+h]y3
% with y1(0)=y1(10)=0, y3(0)=0, y4(0)=1. 
%NOTE: y3= dy1/ dy2(0) and y4=dy2/dy2(0) where d here means partial
%derivative.

%-------------------------------------------------------------------------
function main
clf 
clear all
%% ------ Initialize parameters and variables -------%%
global L hmax h0 h1 h2 h3 h4 h N y20 lam1L ds dxstep errtol dxlam
L = 10;                    % length of habitat
hmax = 1;                  % maximum harvesting
% constructing the control form
h0 = hmax;                 % h(x) = hmax on (s0,s1) where s0=0
h1 = 0;                    % h(x) = 0    on (s1,s2)     
h2 = 0.5;                  % h(x) = 1/2  on (s2,s3) <- the singular arc
h3 = 0;                    % h(x) = 0    on (s3,s4) 
h4 = hmax;                 % h(x) = hmax on (s4,s5) where s5=L=10
h = [h0; h1; h2; h3; h4];  
N = length(h);            
y20 =0.191926627545866;    % initial guess for value y2(0)
y20save = y20;
lam1L = 1.289577736072964;  
ds=L/500;                  % delta =0.02
dxstep =1.e-3;
errtol =1.e-10;
dxlam = dxstep;
% initializing starting guess for each switch
s1 = 1.78;                 % initial guess for switch point s1              
s2 = 2.92;                 % initial guess for switch point s2
s3 = 7.12;                 % initial guess for switch point s3 
s4 = 8.25;                 % initial guess for switch point s4 
s = [s1;s2;s3;s4];
ssave=[s1;s2;s3;s4];       %using for a while loop
olds = s;                  %using for a while loop

%% ------------------bounds on s----------------
%Lower and upper bound for s
format long
% Bounds over each switch. 
% lo and hi are used to indicate that s_i lies in the closed interval 
% [s_i-ds, s_i+ds]
lo = s-ds*ones(length(s),1);
hi = s+ds*ones(length(s),1);
% Each switch is updated through PASA and so is the constraints lo and hi.
oldlo = lo;                 %using for a while loop
oldhi= hi;                  %using for a while loop
%% parameters for pasa and while loop
step = 1;
maxstep = 20;               
%--------- creating bounds of the form As<=bu---------
% This will correspond to constraints being 
% %s_i<=s_{i+1}
n = length(s)-1;
e = ones(n,1);
A = spdiags([e,-e],0:1,n,n+1); 
bl = -inf*ones(n,1);            % bl <= As but bl is unlimitted
bu = zeros(n,1);
%% ================ initializing PASA parameters ===============
pasadata.pasa.PrintStat=1;
pasadata.pasa.grad_tol=errtol;
pasadata.A = A;
pasadata.bl = bl;
pasadata.bu = bu;
pasadata.value = @objective;
pasadata.grad = @grad;
% We run PASA to update the switches s^{k+1}_i and we update the box 
%constraints an iterative number of times until s^{k+1}_i is in the 
%interior of [s^{k}_i-ds,s^{k}_i+ds] for all i=1,2,3,4 or until a maximum 
%number of steps is exceeded. 
%% =========== while loop that is used for solving the problem ==========
while (step<maxstep)
    pasadata.x = olds;
    pasadata.lo = oldlo; 
    pasadata.hi = oldhi;
    [snew, PASAstat] = pasa(pasadata);
    olds= ssave;
    test1 = min(abs(snew'-oldlo));
    test2 = min(abs(snew'-oldhi));
    if min(test1,test2)>1.e-10
        sopt = snew';
        step =maxstep+2;
        break
    else
        ssave = snew';
        olds= snew';
        disp('snew')
        disp(snew)
        [y20,X,Y,yS] = y20update(olds);
        [lam1, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(olds);
        oldlo = olds-ds*ones(length(olds),1);
        disp('oldlo')
        disp(oldlo)
        oldhi = olds+ds*ones(length(olds),1);
        disp('oldhi')
        disp(oldhi)
        step = step+1;
    end
end
%Conditional statement indicating if the maxstep is too small
if step == maxstep
    sopt ='increase maxstep';
end
format long
diary spa4attempt1.log
diary on
disp('sopt')
disp(sopt)
[y20up,X,Y,yS] = y20update(sopt);
disp('y20')
disp(y20)
[lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(sopt);
disp('lam1L')
disp(lam1L)
cost = objective(sopt);
disp('cost')
disp(cost)
diary off


%% testing gradient
% c0 = objective(s) 
% g0 = grad(s)
%     for j = 1:length(s)
%         dt = 0.1 ;
%         for i = 1:8
%             sj = s(j) ;
%             s(j) = s(j) + dt ;
%             y20 = y20save;
%             c1 = objective(s); 
% %             if i ==1 && j==1
% %                 disp(c1)
% %                 disp(y20)
% %                 disp(lam1L)
% %             end
%             g1 = (c1-c0)/dt ;
%             e = g1 - g0(j) ;
%             fprintf (1, 'dt: %10.2e e: %10.2e g1: %25.15e g0: %25.15e\n', ...
%                          dt, e, g1, g0(j)) ;
%             s(j) = sj ;
%             dt = dt/10 ;
%         end
%     end
% 


%% Functions needed for the file
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
% We will use a for loop and cell array commands to save ode45 solution of
% dify over the intervals [S(k),S(k+1)] where S = {0,s1,...,s_N-1,L}
function [X,Y,yS, y1L, y3L] = ysolve(Y20, s)
    S = [0; s; L];
    k = 1;                      %loop counter
    yS = zeros(length(S),5);    %initializing yS- evaluates y at points in S
    yS(1,:) = [0,Y20,0,1,0];    %[y1(0),y2(0),y3(0),y4(0)] 
    intsize = zeros(N,1);       % numb. of mesh pts for solving ode over each interval
    for k = 1 : N
        options = odeset('RelTol',1.e-10,'AbsTol',errtol,'MaxStep', 1.e-2);
        ODE = @(x,y)dify(x,y,h(k));
        [X{k},Y{k}] = ode45(ODE,[S(k),S(k+1)], yS(k,:), options);
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
    while (k<= maxits)
        [X,Y, oldyS, oldy1L, oldy3L] = ysolve(oldy20,s);
        newton = -oldy1L/oldy3L;
        if ( abs(newton) >= 2.e-4 )
            newton = 2.e-4*newton/abs(newton) ;
        end
        newy20 = oldy20+newton ;
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
% This function will call all of the preceeding functions to get a y1
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
% interpolating y1
function lam2int = lam2interpolate(Xlam,Lam)
    k=1; 
    for k = 1:N
        Lamnew = Lam{k};
        lam2{k} = Lamnew(:,2);
        lam2int{k} = spline(Xlam{k},lam2{k});
    end
end
% now we need to getlam2 by calling all of these functions
function lam2int = getlam2(s)
    y1int = gety1(s);
    [lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(s);
    [Xlam, Lam, lam20, lamS] = lamsol2(y1int, s, dxlam);
    [Xlam, Lam, lam20, lamS] = lamsol2(y1int, s, dxlam/2);
    lam2int = lam2interpolate(Xlam,Lam);
end


% %% Functions for Plotting
% %-------- extracting state variable
% function [x, y1plot] = extractstate(X,Y)
%     i = 1; 
%     for i = 1:N
%         ynew = Y{i};
%         Y1{i} = ynew(:,1);
%     end
%     y1plot = vertcat(Y1{:});
%     x = vertcat(X{:});
% end
% %------Extracting lam1 and lam2 for plotting------%
% function [xlam, lambda1, lambda2] = extractcostate(Xlam,Lam)
%     k = N;
%     for k = N:-1:1
%         newlam = Lam{k};
%         Lam1{k} = flip(newlam(:,1)); %flip(lam) = lam values from [s(k), s(k+1)]
%         Lam2{k} = flip(newlam(:,2));
%         xvar{k} = flip(Xlam{k});     %flip(x) = x values from [s(k),s(k+1)]
%     end
%     xlam = vertcat(xvar{:});         % concatenates cell arrays into a vector
%     lambda1 = vertcat(Lam1{:});      % concatenates cell arrays into a vector   
%     lambda2 = vertcat(Lam2{:});      % concatenates cell arrays into a vector
% end
end