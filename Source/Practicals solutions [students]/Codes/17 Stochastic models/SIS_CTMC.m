%%
clear all; clc; clf;

%% Set the Constants 
TOTAL_POPULATION=100; 
TRANSMISSIBILITY=1;   
BIRTH_RATE=0.25; 
RECOVERY_RATE=0.25; 
INTEREVENT_TIME=0.01;
TOTAL_TIME=50; 
NUMBER_OF_SIMUL=10000; 
ZERO_INFECTIVES=2;      

%% Drawing constants
X_AXIS = 50; I_AXIS = 80;

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
b = @(i,dt) beta*i*(N-i)*dt/N;      % Set function handle for simplicity of coding
d = @(i,dt) (bir+gam)*i*dt;
s = @(i,dt) -(b(i,dt)+d(i,dt));

TrMat = zeros(N+1,N+1);              % Set the Transition probability (i=0~N)
TrMat(1,1) = 1; % Zero state is absorbing point.
for i = 2:N+1
  TrMat(i-1,i) = d(i-1,1);
  TrMat(i,i) = s(i-1,1);
  TrMat(i+1,i) = b(i-1,1);
end
TrMat(N+1,N+1) = 1-d(N,1);
Q = TrMat(2:N+1,2:N+1);
TrMat(N+2,:) = [];

%% Simulation
Allsim=cell(NUMBER_OF_SIMUL,1); % Set the simulation 
AllT=cell(NUMBER_OF_SIMUL,1);
last_time_I=zeros(NUMBER_OF_SIMUL,1);
for j = 1:NUMBER_OF_SIMUL
    i=2; 
    clear T; clear simul;
    simul(1) = I0;
    T(1) = -log(rand)/(b(I0,1)+d(I0,1));            % For interevent time T, there is no dt term becase it is deleted
    while(sum(T)<TOTAL_TIME)
        prei=simul(i-1);
        a1 = b(prei,T(i-1))/(b(prei,T(i-1))+d(prei,T(i-1))); 
        a2 = a1 + d(prei,T(i-1))/(b(prei,T(i-1))+d(prei,T(i-1)));
        c = rand;
        if b(prei,T(i-1))+d(prei,T(i-1)) == 0
            simul(i) = simul(i-1);
        else
            if(c<=a1)
                simul(i) = simul(i-1)+1;
            elseif(c>a1 && c<=a2)
                simul(i) = simul(i-1)-1;
            end
        end
        if b(simul(i),1)+d(simul(i),1) == 0
            T(i)=0.01;
        else
            T(i)=-log(rand)/(b(simul(i),1)+d(simul(i),1));
        end
        i=i+1;
    end
    T=cumsum(T);
    T(length(T))=[];
    t=[0 T];
    last_time_I(j)=simul(:,end);
    Allsim{j,1} = simul;
    AllT{j,1} = t;
end

%% Quasi-state stationary points
distribution_of_I=zeros(N+1,1);
for i=0:N
    distribution_of_I(i+1)=sum(last_time_I==i);
end
quasi_numerical=distribution_of_I(2:end)/(NUMBER_OF_SIMUL-distribution_of_I(1));

quasi_eq=@(q) Q*q+(bir+gam)*q(1)*q;
quasi=fsolve(quasi_eq, ones(N,1)/N);

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
figure(2);
hold on;                          % Figure 3.7
for i=1:NUMBER_OF_SIMUL
    t=AllT{i,1}; simul=Allsim{i,1};
    plot(t,simul);
end
plot(td,yd(:,2),'k-','LineWidth',4);
axis([0 X_AXIS 0 I_AXIS]);
xlabel('Time');
ylabel('Number of Infectives');
hold off;

%% Histogram of paths
figure(3)
hold on;
cumDat=[]; T=[];
for i =1: NUMBER_OF_SIMUL
    T = [T AllT{i,1}];
    cumDat = [cumDat Allsim{i,1}];
end
h=histogram2(T, cumDat, 'BinWidth',[TIME_UNIT,INFECTIVE_UNIT],'FaceColor','flat','DisplayStyle','bar3','ShowEmptyBins','on');
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
Y(1)=[];
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