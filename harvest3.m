%% ---------------Problem Description--------------------%
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

%This file uses total variation regularization to obtain  
%a regularized solution to this problem. 
%The regularization tuning parameter is p=10^{-3}

%The numerical scheme used here is detailed in Chapter 3 of
% Regulariztion of Singular Control Problems that Arise in Mathematical
% Biology
% Author: Summer Atkins
% year: 2021
% PhD Thesis, Department of Mathematics, University of Florida
% url https://ufdc.ufl.edu/UFE0057379/00001/pdf

% You will need to download SuiteOPT software in order to run this file
% This software can be accessed on
% https://people.clas.ufl.edu/hager/software/ 
% In our experiments we used the SuiteOPT version 2.0.3. You can use the 
% same link but click software archive to access this older version of the
% software. 
function main
clear all
clf
%%-----------------Begining Script-----------------------%
global L N u0 uL p tol hmax lamV1
L = 10;             % spatial boundary
N = 500;      
u0 =0;              %boundary conditions on state variable u
uL = 0;
p =1.e-3;           % total variation penalty parameter 
tol = 1.e-8;        % KKT tolerance for pasa
hmax = 1;           % upperbound on control
lamV1 = 0;          % transversality condition on lambdaV
%%=======================Parmeters for PASA=======================
%---------defining control variable h----------------------
h = zeros(N,1);
%zeta and iota will be used for calculating the total variation of h 
zeta = zeros(N-1,1); 
iota = zeros(N-1,1);
H = [h; zeta; iota];
%-------------Initialize Matrix A and bounds for PASA -------------%
%we want bl<= AH <= bu 
%A is used in ensuring the decomposition of the total variation
%regularization term. 
% let h(i+1) -h(i) = zeta(i)-iota(i)
%If h(i+1)-h(i)> 0 , then zeta(i) = h(i+1)-h(i) and iota(i) = 0;
%If h(i+1)-h(i)<= 0, then zeta(i) = 0 and iota(i) =-(h(i+1)-h(i));
%We want to find a matrix A and bounds 'bl' and 'bu' that correspond 
% to the above conditional statement. 
%We use spdiags to find A
B1 = [-ones(N-1,1), ones(N-1,1), -ones(N-1,1), ones(N-1,1)];
d1 = [0,1,N, 2*N-1];
A1 = spdiags(B1,d1,N-1,N+2*(N-1));
%A is a  N-1 by 3N-2 matrix that reflects the equations:
%h_{k+1} -h_{k}-zeta_{k}+iota_{k}=0 for all k =1...N-1
%------------bounds bl and bu-------
% bl and bu will be all zeros N-1 vector
bl = zeros(N-1,1);
bu = zeros(N-1,1);
%---------Bounds on Control H: lowU and upU---------%
%we want zeta and iota to be non-negative: 0<= zeta_k, iota_k <= inf
%we want 0<= h_k <= hmax;
lowU = zeros(N+2*(N-1),1);
upU = hmax*ones(N+2*(N-1),1);
upU(N+1:3*N-2)= Inf; %zeta and iota are not bounded from above

%% --------------User defined parameter values for pasa---------------%%
% See pasa_default inside Source/pasa.c for a description of the parameters
%see pproj_default inside PPROJ/Source/pproj.c for a description of the 
% parameters. 
pasadata.pasa.PrintStat=1;      % print statistics for used routines
pasadata.pasa.grad_tol= tol;    % stopping tolerence for PASA
%-----------Setup pasadata-------------%
%Since pasadata.lambda, the mulitpliers associated with the constraint 
%bl <= AU<= bu, will not be given, pasa will set pasadata.lambda as being a 
% zero vector by default.

%% --------------------- Set options for pasa -------------------- %
pasadata.x  = H;                % control
pasadata.lo  = lowU;            % lower bound of the control  
pasadata.hi   = upU;            % upper bound of the control
pasadata.A  = A1;               % sparse matrix used for constraints 
pasadata.bl  = bl;              % lower bound for A1H
pasadata.bu = bu;               % upper bound for A1H
pasadata.value = @obj;          % objective function handle
pasadata.grad = @grad;          % gradient function handle
%% ======================CALLING PASA===================================%
tic
[Hp,PASAstat]=pasa(pasadata);
runtime = toc; 
fprintf('N %g\np %e\ntol %g\nL %g \n', N, p, tol, L)
fprintf('runtime:%e\n', runtime);
%hval=Hp(250);
format long
%disp('Approximate Value of h along what should be singular region');
%disp(hval);
znew = updatestate(Hp);
disp('y2(0)')
disp(znew(1));
lambdanew= updatecostate(znew,Hp);
disp('lambda_1(L)')
disp(lambdanew(2*N-2));
%% approximating switches
disp('Finding switches to regularized control')
k= 1;
%1st switch (h switches from hmax to 0
while(Hp(k)>tol)                            %while h(x) is greater than 0
    k =k+1;
end
% fprintf('k: %e',k);
x = linspace(0,L,N+1);                      % partitioning x 
pswitch1=x(k);                              % computing value of s1
s1=x(k);                                    % switch s1
fprintf('first switch: %e \n',pswitch1);
disp('')
%2nd switch (h switches from 0 to 1/2 (actually approximately 0.5134)
while(Hp(k)<0.5)                            
    k = k+1;
end
% fprintf('k: %e', k);
pswitch2= x(k); 
s2 = x(k);                                  % switch s2 
fprintf('second switch: %e \n', pswitch2);
disp('')
% 3rd switch: h switches from 1/2 to 0
while (Hp(k) > tol)
    k = k+1;
end

pswitch3= x(k);
s3 = x(k);                                 % switch s3

fprintf('third switch: %e \n', pswitch3);
disp('')
%4th switch H switched from 0 to hmax
while (Hp(k)<=tol)
    k = k+1;
end
% fprintf('k: %e', k);
pswitch4 = x(k);
s4 = x(k);                                  % switch s4
fprintf('fourth switch: %e \n', pswitch4);
disp('')
z = updatestate(Hp);
[up,vp]=stateunshuffle(z);
lambdap = updatecostate(z,Hp);
[lamUp,lamVp]= costateunshuffle(lambdap);
lamVp= [lamVp(1);lamVp];
lamUp = [lamUp(1);lamUp];
H1 = psi(Hp);

Jp = obj(Hp);

fprintf('approximated regularized optimal value: %e\n',Jp);
disp('')
p = 0;
format long
Jnp = obj(Hp);

fprintf('approximated unregularized optimal value: %e\n',Jnp);
% Jnp

%% ------Plotting solution for unpenalized problem-------
hp = [Hp(1);Hp(1:N)'];
x = linspace(0,L,N+1);

figure(1)
plot(x,up,'-b','LineWidth', 2);
hold on
plot(x,hp, 'r', 'LineWidth', 2);
ax = gca;
ax.FontSize = 25;
%savefig('harvest4state4.fig')

figure(2) 
plot(x, H1, 'k', 'LineWidth', 2);
hold on
line([s1,s1],[-0.05,0.01],'Color', 'red', 'LineStyle', '--')
hold on 
line([s2,s2],[-0.05,0.01],'Color', 'red', 'LineStyle', '--')
hold on 
line([s3,s3],[-0.05,0.01],'Color', 'red', 'LineStyle', '--')
hold on 
line([s4,s4],[-0.05,0.01],'Color', 'red', 'LineStyle', '--')
hold on
plot(x, zeros(size(x)),'--g','LineWidth',1)
ylim([-0.05,0.01])
xlabel('Location in Habitat','FontSize',20)
ylabel('Switching function at Location x','FontSize',20)
title('Switching function with Penalty p =10^{-3}','FontSize',20)
%savefig('harvestswitch3v2.fig');

%% ---------------------SHUFFLING FUNCTION HANDLES-----------------------%%
%we will need to shuffle and unshuffle state and costate variables
%this is due to our formulization of our discretization
%-------------------------SHUFFLE FUNCTION------------------------------%
%this function will be used for other function handles
    function z = shuffle(x,y) 
        z = [x(:).'; y(:).'];
        z = z(:).';
        z =z';
        %z=[x0,y0,x1,y1,...,xN,yN]
    end
%-------------------------UNSHUFFLE FUNCTION-----------------------------%
%this funtion will be used for other function handles
    function [x,y] =unshuffle(z)
        %input z=[x0,y0,x1,y1,...,xN,yN]
        %output x = [x0,x1,...,xN] and y=[y0,y1,..., yN]
        x = z(1:2:end);
        y = z(2:2:end);
    end
%+++++++++++++++++++SHUFFLING STATE VARIABLES++++++++++++++++++++++++++++%
%will be using this function handle for the function handle involving 
%updating state variables
     function z = stateshuffle(u,v)
        %input: u=[u0,u1,...,uN]' and v=[v0,v1,...,vN]'
        %Output: z=[v0,u1,v1,u2,...,uN-1,vN-1,vN]'
        z = shuffle(u,v);
        z(1)=[];
        z(2*N)=[];
    end
%++++++++++++++++++++UNSHUFFLING STATE VARIABLES++++++++++++++++++++++++++%
%This function will be used a lot in regards to function handles used fro
%updating STATE variables as well as function handles used for computing
%the objective function and for finding the gradient of the objective
%function.
    function [u,v]=stateunshuffle(z)
        %input: z =[v0,u1,v1,u2,...,uN-1,vN-1,vN]' 
        %output: u = [u0,u1,...uN]' and v= [v0,v1,...,vN]'
        newz=[z(1:2*N-1);uL;z(2*N)];
        newz=[u0;newz];
        [u,v] = unshuffle(newz);
    end
%++++++++++++++++++SHUFFLING COSTATE VARAIBLES++++++++++++++++++++++++++++%
 %this function will be used in updating our COSTATE variables
    function lambda = costateshuffle(lamU,lamV)
        %inputs:lamU =[lamU1,lamU2,...,lamUN]' & lamV =[lamV1,...,lamVN]'
        %output: lambda =[lamU1,lamU2,lamV2,lamU3,...,lamVN-1,lamUN,lamVN]'
        %notice that lamV1 is omitted from output lambda
        lambda = shuffle(lamU,lamV);
        lambda(2) =[];
        %If lamU and lamV are N vectors 
        %then output lambda should be a 2N-1 vector
    end
%+++++++++++++++++UNSHUFFLEING COSTATE VARIABLES++++++++++++++++++++++++++%
%this function will be used in updating our COSTATE variables
    function [lamU,lamV] =costateunshuffle(lambda)
        %input: lambda =[lamU1,lamU2,lamV2,lamU3,...,lamVN-1,lamUN,lamVN]'
        %outputs: lamU =[lamU1,lamU2,...,lamUN]' & lamV =[lamV1,...,lamVN]'
        Lambda=[lambda(1);lamV1;lambda(2:2*N-1)];
        [lamU,lamV]=unshuffle(Lambda);
    end
%% ==================FUNCTIONS USED FOR UPDATING STATE=====================%
%For updating state equation u we converted the 2nd order eqn into
%first order system         u' = v;                  
%                           v' = hu-u(1-u);
%with boundary conditoins u(0) = u(L) =0;
%we used an implicit euler scheme for discretizing state variables: 
%u =(u0,u1,...,uN)' where u0=u(0) and uN=u(L)
%v = (v0,v1,...,uN)' 
%control: h = (h1,h2,...,hN)
%for all k=0,1,...N-1
%u{k+1}-uk = dx*(v{k+1})
%v{k+1}-vK = dx*(h{k+1}u{k+1}-u{k+1}+u{k+1}*u{k+1})
%We used implicit because there were no initial conditions held on v 
%in addition through discretization via implicit euler we had
%vN-1 -vN = dx(hN*uN-uN+uN*uN)= 0 since uN= u(L). 
%Even if we used an explicit euler scheme we would have had to use newtons
%method for solving this scheme
% for all k we move terms to left hand side 
%u{k+1}-uk - dx*(v{k+1})=0
%v{k+1}-vK - dx*(h{k+1}u{k+1}-u{k+1}+u{k+1}*u{k+1})=0
%WE'LL USED NEWTON'S MATHOD FOR SYSTEMS
%-------------functions needed for Newton's method------------------------%
    function F = f(z,H)
        % function that evaluates the 2N system w.r.t. states & control
        %inputs:state z=[v0,u1,v1,u2,...,uN-1,vN-1,vN]' and 
        %       control H=[h1,...,hN,zeta1,...,zetaN-1,iota1,...,iotaN-1]';
        %output 2N vector where
        %odd entries:  u_{k+1}-u_{k}-dx*v_{k+1}
        %even entries: v_{k+1}-v_{k}-dx*(h_{k+1}*u_{k+1}-u_{k+1}+u_{k+1}^2 
        dx = L/N; %stepsize
        h = H(1:N);
        [u,v]=stateunshuffle(z);
        for i = 1:N
            fu(i) = u(i+1)-u(i)-dx*v(i+1);
            fv(i) = v(i+1) -v(i) -dx*(h(i)*u(i+1)-u(i+1)+u(i+1)^2);
        end
        F =shuffle(fu,fv);
    end
%---------Computing the Jacobian of F
    function dF = df(z,H)
    %inputs: state z and control H
    %output: 2N by 2N matrix that represents the Jacobian of F
    %will be using spdiags(B,d,m,n) to compute the jacobian
        dx = L/N; %stepsize
        h =H(1:N);
        [u,v] = stateunshuffle(z);
        b1 = -1.*ones(2*N,1); %first diagonal
        b3 = [ones(2*N-1,1);-dx]; %third diagonal
        for i = 1:N-1
            s(i) = -dx*(h(i)-1+2*u(i+1));
        end
        s(N) =1;
%        s
        b2= shuffle(zeros(N,1),s); %second diagonal (true diagonal)
        b4= shuffle(-dx.*ones(N,1),zeros(N,1)); %fourth diagonal
        B = [b1,b2,b3,b4];
         d= [-1;0;1;2];
        dF= spdiags(B,d,2*N,2*N);
        %full(dF)
    end
%-----------------------UPDATING STATE EQUATIONS-------------------------%
    function zstar = updatestate(H)
        dx = L/N; %stepsize
        %define variables u and v
        u = 2.*ones(N+1,1);
        u(1) = u0;
        u(N+1)=uL;
        v= ones(N+1,1);
        oldz = stateshuffle(u,v);
        %variables for for newton's method
        errtol = tol; %you can set it to whatever
        maxits = 100; 
        k = 1;
        while (k <= maxits)
                F= f(oldz,H);
                dF= df(oldz,H);
                y = -dF\F;
                newz = oldz+y;
                if norm(y)< errtol
                   zstar= newz;
                   % fprintf('terminated at k=: %e', k);
                   break
                end
                k= k+1;
                oldz = newz;
        end
        if k == maxits+1
           zstar= 'Maximum number of iterations exceeded';
           disp(zstar);
        end
    end
%% =======================UPDATING COSTATE EQUATIONS=====================%
    %Euler's  implicit scheme to adjoint eqns gave us a linear system
    %So the function handle updatecostate() will be solve the linear system
    %to create the matrix that represents the linear system-use spdiags
    function lambdastar= updatecostate(z,H);
        dx = L/N; % stepsize
        h = H(1:N);
        %define lambda
        lamU= ones(N,1);
        lamV= 2.*ones(N,1);
        lambda = costateshuffle(lamU,lamV);
        [lamU,lamV] = costateunshuffle(lambda);
        z = updatestate(H);
        [u,v]=stateunshuffle(z);
        %using spdiags to create matrix for linear system
        %preparing to define matrix B that is used for spdiags
       for i =1:N
            r(i) = dx*(h(i)-1+2*u(i+1));
        end
        r = r';
        x = r;
        x(1) = -1;%first entry of b3 is -1 
        x(N) = -1;%last entry of b4 is -1
        b1 = shuffle(zeros(N,1),dx.*ones(N,1));
        b1(2*N) = [];
        b2 = [dx;-1.*ones(2*N-4,1);dx;dx];
        b3 = shuffle(x, zeros(N,1));
        b3(2*N) = [];
        b4 = ones(2*N-1,1);
        B = [b1,b2,b3,b4];
        d = [-2,-1,0,1];
        A = spdiags(B,d,2*N-1,2*N-1);
        val(1) = dx*h(1)-lamV(1)*r(1);
        for j =2:N-1
            val(j) = dx*h(j);
        end
        val(N) =0;
        val =val';
        val =shuffle(val,zeros(N,1));
        val(2) = -lamV(1);
        val(2*N) = [];
        lambdastar = A\val;
    end
%% ------------------------- objecitve function handle--------------------%
     function J = obj(H)
            dx = L/N; %stepsize
            h = H(1:N);
            zeta = H(N+1:2*N-1);
            iota = H(2*N:3*N-2);
            %update state
            z = updatestate(H);
            [u,v] = stateunshuffle(z);
            lambda = updatecostate(z,H);
            int = 0 ;
            for i = 1:N
                int = int-dx*h(i)*u(i+1);
            end
            totalvar= sum(zeta+iota);
            J = int+ p*totalvar;
     end
 %% ----------------------gradient of objective function handle-----------
    function g = grad(H) 
        dx = L/N; %stepsize
        h = H(1:N);
        zeta = H(N+1:2*N-1);
        iota = H(2*N:3*N-2);
        %update state
        z = updatestate(H);
       [u,v] = stateunshuffle(z);
       %update costate
       lambda = updatecostate(z,H);
       [lamU,lamV] =costateunshuffle(lambda);
       for i = 1: N
           g(i) =dx*u(i+1)*(lamV(i)-1);
       end
       g = g';
       pen = p*ones(2*N-2,1);
       g = [g;pen];
    end
%% ======================switching function =========================
    function H1 = psi(H)
        dx = L/N;
         z = updatestate(H);
       [u,v] = stateunshuffle(z);
       %update costate
       lambda = updatecostate(z,H);
       [lamU,lamV] =costateunshuffle(lambda);
       lamV = [lamV(1);lamV];
       for i = 1:N+1
       H1(i)=(lamV(i)-1)*u(i);
       end
       H1 = H1';
    end
        
end