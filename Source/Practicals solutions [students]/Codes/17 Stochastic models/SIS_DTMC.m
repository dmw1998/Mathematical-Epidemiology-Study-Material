%%
clear; clc; clf;

%% Set the Constants 
TOTAL_POPULATION=100; 
TRANSMISSIBILITY=0.25;   
BIRTH_RATE=0.25; 
RECOVERY_RATE=0.25; 
INTEREVENT_TIME=0.01;
TOTAL_TIME=50; 
NUMBER_OF_SIMUL=5; 
ZERO_INFECTIVES=20;      

%% Drawing constants
X_AXIS = 50; I_AXIS = 30;

%% For simplicity of coding
TIME_UNIT = 0.5; INFECTIVE_UNIT=1; % For histogram

N = TOTAL_POPULATION; 
beta = TRANSMISSIBILITY; 
bir = BIRTH_RATE; 
gam = RECOVERY_RATE; 
dt = INTEREVENT_TIME;
t=0:dt:TOTAL_TIME; 
I0 = ZERO_INFECTIVES;                    

%% Stochastic SIS
b = @(i) beta*i*(N-i)*dt/N;      % Set function handle for simplicity of coding
d = @(i) (bir+gam)*i*dt;
s = @(i) 1-b(i)-d(i);

T = zeros(N+1,N+1);              % Set the Transition probability (i=0~N)
T(1,1) = 1; % Zero state is absorbing point.
for i = 2:N+1
  T(i-1,i) = d(i-1);
  T(i,i) = s(i-1);
  T(i+1,i) = b(i-1);
end
T(N+1,N+1) = 1-d(N);
T(N+2,:) = [];

%% Simulation
simul = zeros(NUMBER_OF_SIMUL,length(t));   % Set the simulation 
simul(:,1) = I0;
for j = 1:NUMBER_OF_SIMUL
  for i = 2:length(t)
     if(simul(j,i-1)==0)
        a1=1; a2=1;
     elseif(simul(j,i-1)==N)
        a1=T(N+1,N+1); a2=1; 
     elseif(simul(j,i-1)~=0 && simul(j,i-1)~=N)
        k=simul(j,i-1)+1;
        a1=T(k,k); a2=T(k-1,k)+T(k,k);
     end
     
     c=rand;
     if(c<=a1)      
       simul(j,i)=simul(j,i-1);
     elseif(c>a1 && c<=a2)
       simul(j,i)=simul(j,i-1)-1;
     elseif(c>a2 && c<=1)
       simul(j,i)=simul(j,i-1)+1;
     end
  end 
end

%% Quasi-state stationary points
last_time_I=simul(:,end);
distribution_of_I=zeros(N+1,1);
for i=0:N
    distribution_of_I(i+1)=sum(last_time_I==i);
end
quasi_numerical=distribution_of_I(2:end)/(NUMBER_OF_SIMUL-distribution_of_I(1));

state_pdf = zeros(N+1,length(t));     % Calculate the probability of each state
state_pdf(I0+1,1) = 1;
for i = 2:length(t)
  state_pdf(:,i) = T*state_pdf(:,i-1);
end

% Note that 'state_pdf' includes 0 states.
quasi=state_pdf(2:end,2:end)./(1-state_pdf(1,2:end));

figure(1);
hold on;
plot(1:N,quasi(:,end),'k-','LineWidth',1);
plot(1:N,quasi_numerical, 'd');
legend(["Theory","Simulation"]);
xlabel("Quasi-stationary state");
ylabel("Probability");
hold off;

%% Deterministic SIS
I(1) = ZERO_INFECTIVES; S(1)=N-I(1);        % Comparing for Deterministic Model
endtime=TOTAL_TIME;                          
divide=500;                            
td = linspace(0,endtime,divide); h=td(2)-td(1);
y = [S(1); I(1)];     
yd = zeros(1,2);       
yd(1,:) = y';   
f = @(y) [-beta*y(2)/N bir+gam;beta*y(2)/N -(bir+gam)]*y;       
for i=2:divide 
  k1 = f(y);          
  k2 = f(y+k1*h/2);   
  k3 = f(y+k2*h/2);     
  k4 = f(y+k3*h);     
  y = y + h*(k1+2*k2+2*k3+k4)/6;  
  yd(i,:) = y';       
end

%% Plot realizations
figure(2)
hold on;                          
for j=1:NUMBER_OF_SIMUL
    plot(t,simul(j,:));
end
plot(td,yd(:,2),'k-','LineWidth',4);
% plot(td,N-yd(:,2),'r--'); % Susceptible..
axis([0 X_AXIS 0 I_AXIS]);
xlabel('Time');
ylabel('Number of Infectives');
hold off;

%% Histogram of paths
figure(3)
hold on;
cumDat=[]; T=[];
for i =1: NUMBER_OF_SIMUL
    T = [T t];
    cumDat = [cumDat simul(i,:)];
end
h=histogram2(T, cumDat, 'BinWidth',[TIME_UNIT,INFECTIVE_UNIT],...
    'FaceColor','flat','DisplayStyle','bar3','ShowEmptyBins','on');
colorbar
plot(td,yd(:,2),'k-','LineWidth',4);
axis([0 X_AXIS 0 I_AXIS]);
xlabel('Time');
ylabel('Number of Infectives');
hold off;

figure(4)
hold on;
t2=h.XBinEdges;
t2(end)=[];
Y=h.YBinEdges;
Y(end)=[];
count = h.Values;
siz = size(count);
maxData = max(count');
for i= 1:siz(1)
    count(i,:)=count(i,:)/maxData(i);
end
b=bar3(Y,count',1);
for k=1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'flat';
end
colorbar
Xdat = get(b,'XData');
for ii=1:length(Xdat)
    Xdat{ii} = Xdat{ii}/(1/TIME_UNIT);
    set(b(ii),'XData',Xdat{ii});
end

axis([0 X_AXIS 0 I_AXIS]);
xlabel('Time Step');
ylabel('Number of Infectives');
hold off;