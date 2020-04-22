clc
clear all;
close all;

data_choose = "data3";
Coor = load('data\' + data_choose + '.txt');
N = length(Coor(:,1));

%the final regression result is affected by initial values
if data_choose == "data4"
    [x0,y0,r] = circle_es(Coor(1,1),Coor(1,2),Coor(30,1),Coor(30,2),Coor(60,1),Coor(60,2));
    epsilon = 0.05;
elseif data_choose == "data3"
    [x0,y0,r] = circle_es(Coor(20,1),Coor(20,2),Coor(50,1),Coor(50,2),Coor(70,1),Coor(70,2));
    epsilon = 0.3;
elseif data_choose == "data2"
    [x0,y0,r] = circle_es(Coor(10,1),Coor(10,2),Coor(40,1),Coor(40,2),Coor(70,1),Coor(70,2));
    epsilon = 0.1;
elseif data_choose == "data1"
    [x0,y0,r] = circle_es(Coor(5,1),Coor(5,2),Coor(50,1),Coor(50,2),Coor(100,1),Coor(100,2));
    epsilon = 0.01;
elseif data_choose == "data0"
    [x0,y0,r] = circle_es(Coor(20,1),Coor(20,2),Coor(50,1),Coor(50,2),Coor(70,1),Coor(70,2));
    epsilon = 0.05;
else
    disp('input error');
end
    

%%    
%%%%%%%%%%%%%%%%%%%%%%
%%%%--f(X)=H*X--%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%%
%the construction of f(X)


X = [x0 x0^2 y0 y0^2 r^2 0.001*ones(1,N) 0.001*ones(1,N)]';
%the variable that you should initialize
%the size of X is (2N+5)*1

C = 10;
%regularization constant

H = C*[zeros(1,5) ones(1,2*N)];
%a coefficient matrix
%the size of H is 1*(2N+5)
%%

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%--A*X<=b--%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%%
%the construction of A
%the size of A is 4N*(2N+5)

Coor_x = Coor(:,1);
%horizontal ordinate of all samples

Coor_y = Coor(:,2);
%longitudinal coordinates of all samples

A1 = [-2*Coor_x ones(N,1) -2*Coor_y ones(N,1) -1*ones(N,1) -1*diag(diag(ones(N,N))) zeros(N,N)];

A2 = [2*Coor_x -1*ones(N,1) 2*Coor_y -1*ones(N,1) ones(N,1) zeros(N,N) -1*diag(diag(ones(N,N)))];

A3 = [zeros(N,5) -1*diag(diag(ones(N,N))) zeros(N,N)];

A4 = [zeros(N,5) zeros(N,N) -1*diag(diag(ones(N,N)))];

A = [A1;A2;A3;A4];

%the construction of b
%the size of b is 4N*1

%epsilon = 0.05;

b1 = epsilon*ones(N,1) - Coor_x.^2 -Coor_y.^2;

b2 = -1*b1;

b3 = zeros(N,1);

b4 = b3;

b = [b1;b2;b3;b4];
%%

%%%%%%%--Augmented Lagrangian--%%%%%%%
%----Lp(X,Y) = f(X) + Y'(AX-b-m) + rho/2*||AX-b-m||^2----%

%Y is the Lagrange multiplier
%the size of Y is 4N*1

%m = AX-b
%m<=0
%the size of m is 4N*1

Y = 5*ones(4*N,1);

rho = 10;
%the penalty factor

% [X,fval_linear] = linprog(H,A,b,[],[]);


%%

%%%%%%%--iterative process--%%%%%%%%

%%%%%--iteration k+1--%%%%%
% m(k+1) = min(AX(k)-b+Y(k)/rho,0), for every component of vector m
% 
% X(k+1) = (A'*A)^(-1)*((rho*A'*(b+m(k+1)) - H' - A'*Y(k))/rho)
% 
% Y(k+1) = Y(k) + rho*(A*X(k+1) - b - m(k+1))

for k = 1:1:1000
    m = zeros(4*N,1);
    m_temporary = A*X-b+Y./rho;
    small = (m_temporary<=zeros(4*N,1));
    min_flag = find(small == 1);
    m(min_flag,:) = m_temporary(min_flag,:);
    
    X = (A'*A)\((rho.*A'*(b + m) - H' - A'*Y)./rho);
    
    Y = Y + rho.*(A*X - b -m);
end


clc
close all;
figure(1);
plot(Coor_x,Coor_y,'r.');
hold on;

radius =sqrt(epsilon + X(5) - X(2) -X(4) + X(1)^2 + X(3)^2);
rectangle('Position',[X(1)-radius,X(3)-radius,2*radius,2*radius],'Curvature',[1 1],'linewidth',1); 
axis equal




