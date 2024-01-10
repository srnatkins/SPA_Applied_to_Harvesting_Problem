 %%---------------------------Intro.------------------------- 
% This script file uses PASA to solve for Neubert's harvesting problem via 
% a switch point algorithm in the case where the habitat length L=10 
% and the parameter hmax=1. In this file we assume the control form has 
% 6 switches with no singular arc.

% When using the switches obtained from h_{6,sing}^* for h_{6,bang}, the 
% resulting state y1 satisfies the boundary conditions (after 
% applying Newton's method), y1 is negative on certain regions. With this
% in mind we need to adjust the switches to where the resulting state is a
% physically relevant solution. 

% We use PASA to update the switches that came from h_{6,sing}^* one time. 
% The box constraints used in updating the switches is [s_i-0.04, s_i+0.04] 
% for i=1,2,...,6. After this update, the resulting state solution is
% non-negative everywhere. 

% You will need to download SuiteOPT software in order to run this file
% This software can be accessed on
% https://people.clas.ufl.edu/hager/software/ 
% In our experiments we used the SuiteOPT version 2.0.3. You can use the 
% same link but click software archive to access this older version of the
% software. 

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

%Note in our file y1(x)=u(x), y2(x)=u'(x). Additionally for our state
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
global L hmax h0 h1 h2 h3 h4 h5 h6 h N y20 lam1L ds dxstep errtol dxlam
L = 10;                    % length of habitat
hmax = 1;                  % maximum harvesting
% constructing the control form
h0 = hmax;                 % h(x) = hmax on (s0,s1) where s0=0               
h1 = 0;                    % h(x) = 0    on (s1,s2)     
h2 = hmax;                 % h(x) = hmax on (s2,s3) 
h3= 0;                     % h(x) = 0    on (s3,s4) <- No singular arc 
h4 = hmax;                 % h(x) = hmax on (s4,s5)
h5 = 0;                    % h(x) = 0    on (s5,s6)
h6 = hmax;                 % h(x) = 0    on (s6,s7) where s7=L=10
h = [h0; h1; h2; h3; h4; h5; h6];         
N = length(h);      
y20 = 0.192080841873656;
y20save = y20;
lam1L =  1.301535317935051; 
ds = L/250;                 %ds=0.04
dxstep = 1.e-3;
errtol =1.e-10;
dxlam = dxstep;

% initializing starting guess for each switch
%These came from using SPA to improve the switches of h_{6,sing}. 
s1=1.768121762005436;  
s2=3.103103235269388;   
s3=3.357266684882776;   
s4=6.642733319818394;
s5=6.896896766723736;
s6=8.231878237927704;


s = [s1;s2;s3;s4;s5;s6];
[y20up,X,Y,yS] = y20update(s);
[x, y1plot1] = extractstate(X,Y);
disp('y20')
disp(y20)
% if one were to plot the state solution they would see that y1 is negative
% for certain regions of the habitat, which is not physically relevant. 
% figure(1)
% plot(x,y1plot1)

%% ------------------bounds on s----------------
%Lower and upper bound for s
format long
% Bounds over each switch. 
%lo and hi are used to indicate that s_i lies in the closed interval [s_i-ds, s_i+ds]
lo = s-ds*ones(length(s),1);
hi = s+ds*ones(length(s),1); 


%--------- creating bounds of the form As<=bu---------
% This will correspond to constraints being 
% %s_i<=s_{i+1}
n = length(s)-1;
e = ones(n,1);
A = spdiags([e,-e],0:1,n,n+1); 
bl = -inf*ones(n,1);            % bl <= As but bl is unlimitted
bu = zeros(n,1);

%% ================ setting PASA parameters ===============
pasadata.pasa.PrintStat=1;
pasadata.pasa.grad_tol=errtol;
pasadata.A = A;
pasadata.bl = bl;
pasadata.bu = bu;
pasadata.value = @objective;
pasadata.grad = @grad;
pasadata.x = s;
pasadata.lo = lo; 
pasadata.hi = hi;
[snew, PASAstat] = pasa(pasadata);
snew=snew';
disp('snew')
disp(snew)
[y20,X,Y,yS] = y20update(snew);
[lam1, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(snew);
format long
[x, y1plot] = extractstate(X,Y);
disp('y20')
disp(y20)
[lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(snew);
disp('lam1L')
disp(lam1L)
figure(2)
plot(x,y1plot)
cost = objective(snew);
disp('cost')
disp(cost)

%% Functions needed for file

%% ------Solving for the states--------
%---------------dify - setting up ode for state equations ------------%
function  dy = dify(x,y,h)
    dy(1) = y(2);
    dy(2) = -y(1)*(1-y(1))+h*y(1);
    dy(3) = y(4);
    dy(4) = (2*y(1)-1+h)*y(3);
    dy(5) = -h*y(1);          %RHS is integrand of cost functional 
    dy = dy';
end
%---------------------------------------------------------------------%
%------------------Using ode45 for solving dify--------------%
% we will use a for loop and cell array commands to save ode45 solution of
% dify over the intervals [S(k),S(k+1)] where S = {0,s1,...,s_N-1,L}
function [X,Y,yS, y1L, y3L] = ysolve(Y20, s)
    S = [0; s; L];
    k = 1;                          %loop counter
    yS = zeros(length(S),5);        %yS- evaluates y at points in S
    yS(1,:) = [0,Y20,0,1,0];        %[y1(0),y2(0),y3(0),y4(0)] 
    for k = 1 : N
        options = odeset('RelTol',1.e-10, 'AbsTol', errtol,'MaxStep', 1.e-2);
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
%++++++++Using Newton's Methods to Satisfy Boundary Condition++++++++++++%
function[y20up,X,Y,yS] = y20update(s)
    tol = 1.e-12;
    maxits = 1000;
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
    if k == maxits +1
        y20up = 'Maximum number of iterations exceeded';
        disp(y20up);
    end
end

%% ===================Interpolating y1 for the costate===================
%---------------------------State Equations----------------------------
% this may be unnecessary but it should be a faster ode to solve in
% comparison to the dy = dify(x,y,h) function since it has less state
% variables
function dy2 = dify2(x,y,h)
    dy2(1) = y(2);
    dy2(2) = -y(1)*(1-y(1))+h*y(1);
    dy2 = dy2';
end
%------------------------------------------------------------------------
%-----------------------Using ODE45 to solve dy2----------------------
function [X,Y,yS,y1L] =ysolve2(Y20,s,dx)
    S = [0; s; L];
    k = 1;                      %loop counter
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
% ---------------------------Interpolating y1--------------------------
% function will use spline commands to interpolate
function y1int = y1interpolate(X,Y)
    k = 1;
    for k = 1:N
        Ynew = Y{k};
        y1{k} = Ynew(:,1);
        y1int{k} = spline(X{k},y1{k});
    end
end
% ---------------------- Get y1 for costate variables --------------------
% This function will call the proceeding functions to where the input is
% just the vector s and the output is the interpolating variable y1
function y1int = gety1(s)
    %Use Newton's method to get good y20
    [y20up,X,Y,yS] = y20update(s);   
    %Integrate the state equations using ode45 with  MaxStep =dxstep
    [X,Y,yS,y1L] = ysolve2(y20up,s,dxstep);     
     %Integrate the state equations using ode45 with  MaxStep =dxstep
    [X,Y,yS,y1L] = ysolve2(y20up,s,dxstep/4);   
    %Interpolate y1 using spline
    y1int = y1interpolate(X,Y);                 
end
%--------------------------------------------------------------------------
%==========================================================================
%% ===================Functions for the Costate==========================
%-------costate - function sets up ode system for costate------- %%
function dlambda = costate(x,lambda, y1, h)
    dlambda(1) = h - lambda(2)*(-1+2*ppval(y1,x)+h);
    dlambda(2) = -lambda(1);
    dlambda(3) = (1-2*ppval(y1,x)-h)*lambda(4);
    dlambda(4) = -lambda(3);   
    dlambda = dlambda'; 
end
%---------------Solving for the ode system using ode45 ------------------%
function [Xlam, Lam, lam20, lam40, lamS] = lamsol(Lam1L,s)
    y1int = gety1(s);       %interpolating y1
    S = [0;s;L];
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
%-----------------Using Newton's Method for Updating lam1L ---------------%
function [lam1Lup, Xlam, Lam, lam20, lam40, lamS] = updatelam1L(s)
    maxits = 5;
    tol = 1.e-12;
    j = 1;
    oldlam1L = lam1L;
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
%% -----Functions for Computing the Cost, Hamiltonian, and Gradient------
% ------------------Computing cost function -----------------
function cost = objective(s)
    [y20up, X, Y,yS] = y20update(s);
    cost = yS(end,end);
end
% ------------------- Computing Hamiltonian -------------------
function H = Ham(y1, y2, lam1, lam2,h)
    lam3 = 1;
    H = lam1*y2+lam2*(-y1*(1-y1))+lam2*h*y1+lam3*(-h*y1);
end
% -----------Computing Gradient of Cost wrt switching points------------
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
% -------------------------------------------------------------------
%% --------------- Interpolating the lam2  ----------------
% --------------ode system
function dlambda2 = costate2(x, lambda, y1,h)
    dlambda2(1) = h - lambda(2)*(-1+2*ppval(y1,x)+h);
    dlambda2(2) = -lambda(1);
    dlambda2 = dlambda2';
end
% solving the ode system
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



% %% Functions for plotting
% % optional functions
% %-----------Extracting x & y1 from cell array for plotting-----------%
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
% 

end