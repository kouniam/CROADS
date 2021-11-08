%------------------------------------------------------%
%--------------- DICE99 - ECONOMIC PART ---------------%
%------------------------------------------------------% 

close all;
clear all;

dt = 1;                % time step (in years) 
yr1 = 1;                % first year
yr2 = 350;              % last year 
timvec = yr1:dt:yr2;    % time vector 
nt = length(timvec);    % nr of time steps (which are years)

%------------------ list of variables -----------------%

% *A(t) % economic efficiency or "total factor productivity"
% *c(t) % consumption per capita
% *C(t) % total consumption C(t)=c(t)L(t)
% *g_A(t) % determines time dependence of economic efficiency A
% *g_pop(t) % rate of population increase
% *I(t) % investment per year
% *K(t) % capital
% *L(t) % global population (proportional to labour input)
% *Q(t) % economic output per year
% *R(t) % discount factor arising from pure time preference
% U(t) % current fertility, depending on per capita consumption c and population number L
% W % social welfare function
% *rho(t) % rate of pure time preference

%------------------- model parameters ------------------%

A0 = 0.01685;       % economic efficiency
g_A0 = 0.0038;      % init. value of g_A
g_pop0 = 0.0157;    % init. pop. growth rate
g_rho = 0.0025719;  % reduced impatience over time
K0 = 47;            % init. capital
L0 = 5632.7;        % init. global population
q = 0.2510;         % init. saving rate
lambda = 0.3;       % determines how economic output depends on capital and labour input
delta_A = 1e-9;     % determines changes in g_A
%delta_K = 0.1;      % depreciation of capital at dt=10

delta_K = 0.0633+0.0012*dt+0.00025*(dt^2); % time dependant delta_K

delta_pop = 0.0222; % slowdown of rate of population change

rho0 = 0.03;        % init. rate of pure time preference rho
%rho0 = 0.001;       % Stern Report rate 

% size = 10/0.5;
% 
% i=1;    %indexation number
% W = length(size);
% x_i = length(size);
%  
% for x=0:0.05:10;

%------------------------ initialization values ------------------------%

g_pop(1) = g_pop0;
rho(1) = rho0;
g_A(1) = g_A0;

cumP(1) = (1/(1+rho(1)))^dt;
R(1) = (1/(1+rho(1)))^dt;

L(1) = L0;
A(1) = A0;
K(1) = K0;

Q(1) = A(1)*(K(1)^lambda)*(L(1)^(1-lambda));

I(1) = q*Q(1);

c(1) = (1-q)*(Q(1)/L(1));

C(1) = c(1)*L(1);

c(1) = C(1)/L(1);

w(1) = L(1)*log(c(1)*(1e6))*R(1);

% %--------------- damage ---------------%

% D = zeros(1,nt);
%  
% D(199) = 30; % $10^12

% %--------------------------------------%

for s=2:nt;
        
    rho(s) = rho(1)*exp(-g_rho*s);          % pure rate of time preference
    cumP(s) = cumP(s-1)*(1/(1+rho(s)))^dt;  % comulative product
        
end

for t=2:nt;
    
    R(t) = cumP(t); % discount factor arising from pure time preference
    
    %------------------------ population growth ------------------------%
      
    g_pop(t) = integral(@(s) g_pop(1)*exp(-delta_pop*s),0,t);
        
    L(t) = L(1)*exp(g_pop(t));

    %------------------- total factor productivity ---------------------%

    g_A(t) = integral(@(s) g_A(1)*exp(-delta_A*s),0,t);
    
    A(t) = A(1)*exp(g_A(t));
    
    %------------------------- capital equation ------------------------%

    K(t) = K(t-1).*((1-delta_K).^dt)+dt.*I(t-1);
    
    Q(t) = A(t).*(K(t).^lambda).*(L(t).^(1-lambda));
    
    I(t) = q*Q(t);
    
    c(t) = (1-q)*(Q(t)/L(t));
    
    C(t) = c(t)*L(t);
    
%     --- corrective damage recalculation ---%
%     
%     C(t) = C(t)-D(t);
%         
%     c(t) = C(t)/L(t);
%         
%     ---------------------------------------%

    w(t) = L(t)*log(c(t)*(1e6))*R(t);
    
end

dQ = diff(Q);

G = dQ./Q(1:349);

W = sum(w)  % total wellfare

% W(i) = sum(w);
% x_i(i) = x;
%  
% i=i+1;
%  
% end

figure(1)

subplot(2,2,1);
plot(timvec,K,'Linewidth',2.5)
title('Capital')
xlabel('time [years since 1995]')
ylabel('Capital [10e12$]')

subplot(2,2,2);
plot(timvec,Q,'Linewidth',2.5)
title('Output')
xlabel('time [years since 1995]')
ylabel('Output [10e12$]')

subplot(2,2,3);
plot(timvec,c,'Linewidth',2.5)
title('Per Capita Consumption')
xlabel('time [years since 1995]')
ylabel('Consumption [10e6$]')

subplot(2,2,4);
plot(timvec,w,'Linewidth',2.5)
title('Discounted Utility')
xlabel('time [years since 1995]')
ylabel('Discounted Utility')

figure(2)
plot(timvec,rho,'Linewidth',2.5)
title('Economic Growth')
xlabel('time [years since 1995]')
ylabel('G')

% figure(2)
% plot(timvec(2:349),G(2:349),'Linewidth',2.5)
% title('Economic Growth')
% xlabel('time [years since 1995]')
% ylabel('G')
% 
% figure(3)
% plot(timvec,C,'Linewidth',2.5)
% title('Total Consumption')
% xlabel('time [years since 1995]')
% ylabel('C')

% Wo = (2.497609051163105e6)*ones(1,401); % maximum wellfare function 1.5.a
% Wo = (2.497611138417322e6)*ones(1,401); % maximum wellfare function 1.5.b
% Wo = (2.958405996383646e7)*ones(1,size+1); % maximum wellfare function 1.5.c
% 
%  
% hold on
% 
% plot(x_i,W,'-')
% plot(x_i,Wo,'-r')
% plot(5.947,2.958405996383646e7,'r*')
% title('Accumulated Welfare function vs X')
% xlabel('X')
% ylabel('W')

% [val, idx] = max(W);
% hold on
% plot(q_i(idx),val,'r*')

% hold on
% plot(timvec,gdpA,'Linewidth',2.5)
% plot(timvec,gdp2A,'r','Linewidth',2.5)
% title('GDP comparison (A vs 2*A)')
% xlabel('time [years since 1995]')
% ylabel('Output [10e12$]')
