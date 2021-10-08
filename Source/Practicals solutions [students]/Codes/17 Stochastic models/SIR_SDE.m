%%
clear all; clc; clf;

%% Set the Constants 
TOTAL_POPULATION=100; 
TRANSMISSIBILITY=1;   
BIRTH_RATE=0; 
RECOVERY_RATE=0.5; 
INTEREVENT_TIME=0.01;
TOTAL_TIME=50; 
NUMBER_OF_SIMUL=1000; 
ZERO_INFECTIVES=2;      
ZERO_SUSCEPTIBLES=TOTAL_POPULATION-2;

%% Drawing constants
X_AXIS = 50; I_AXIS = 50;

%% For simplicity of coding
TIME_UNIT = 0.5; INFECTIVE_UNIT=1; % For histogram

N = TOTAL_POPULATION; 
beta = TRANSMISSIBILITY; 
bir = BIRTH_RATE; 
gam = RECOVERY_RATE; 
dt = INTEREVENT_TIME;
t=0:dt:TOTAL_TIME; 
I0 = ZERO_INFECTIVES;      
S0 = ZERO_SUSCEPTIBLES;

%% Stochastic SIR
f1 = @(s,i) beta*i*s*dt/N;      % Set function handle for simplicity of coding
f2 = @(s,i) gam*i*dt;
f3 = @(s,i) bir*i*dt;
f4 = @(s,i) bir*(N-s-i)*dt;
f5 = @(s,i) 1-(f1(s,i)+f2(s,i)+f3(s,i)+f4(s,i));
E = @(s,i) [(-1)*f1(s,i)+(1)*f3(s,i)+(1)*f4(s,i); (1)*f1(s,i)+(-1)*f2(s,i)+(-1)*f3(s,i)];
V = @(s,i) [(-1)^2*f1(s,i)+1^2*f3(s,i)+1^2*f4(s,i) (-1)*1*f1(s,i)+1*(-1)*f3(s,i); ...
    (-1)*1*f1(s,i)+1*(-1)*f3(s,i) 1^2*f1(s,i)+(-1)^2*f2(s,i)+(-1)^2+f3(s,i)];

Data = zeros(NUMBER_OF_SIMUL,length(t));      % Calculate the sample path by SDE of SIR Model
for j=1:NUMBER_OF_SIMUL
    X=[S0; I0];
    Data(j,1) = I0;
    for i=2:length(t)
        [B, D]=eig(V(X(1),X(2)));
        sqrtV=B*sqrt(D)*B^-1;
        normal = [normrnd(0,sqrt(dt)); normrnd(0,sqrt(dt))];
        X=X+E(X(1),X(2))+sqrtV*normal;
        if(X(2)<=0 || Data(j,i-1)==0)
            X(2)=0;
        end
        Data(j,i) = X(2);
    end
end

%% Quasi-stationary distribution (numerical)
last_time_I=Data(:,end);
distribution_of_I=zeros(N+1,1);
for i=0:N
    distribution_of_I(i+1)=sum(last_time_I>=i & last_time_I<i+1);
end
quasi=distribution_of_I(2:end)/(NUMBER_OF_SIMUL-distribution_of_I(1));

figure(1);
plot(1:N,quasi,'kx:','LineWidth',1);
xlabel("Quasi-stationary state");
ylabel("Probability");

%% Deterministic SIR
I(1) = ZERO_INFECTIVES; S(1)=N-I(1);        % Comparing for Deterministic Model
endtime=TOTAL_TIME;                          
divide=500;                            
td = linspace(0,endtime,divide); h=td(2)-td(1);
y = [S(1); I(1)];     
yd = zeros(1,2);       
yd(1,:) = y';   
f = @(y) [bir*N;0] + [-beta*y(2)/N-bir 0;beta*y(2)/N -(bir+gam)]*y;       
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
hold on;                       % Figure 3,10
for j=1:NUMBER_OF_SIMUL
    plot(t,Data(j,:));
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
    T = [T t];
    cumDat = [cumDat Data(i,:)];
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

plot(td,yd(:,2),'k-','LineWidth',4);
axis([0 X_AXIS 0 I_AXIS]);
xlabel('Time Step');
ylabel('Number of Infectives');
hold off;