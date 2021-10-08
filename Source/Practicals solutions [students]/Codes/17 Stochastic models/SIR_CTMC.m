%%
clear all; clc; clf;

%% Set the Constants 
TOTAL_POPULATION=100; 
TRANSMISSIBILITY=1;   
BIRTH_RATE=0.25; 
RECOVERY_RATE=0.25; 
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
f1 = @(s,i,dt) beta*i*s*dt/N;      % Set function handle for simplicity of coding
f2 = @(s,i,dt) gam*i*dt;
f3 = @(s,i,dt) bir*i*dt;
f4 = @(s,i,dt) bir*(N-s-i)*dt;
f5 = @(s,i,dt) 1-(f1(s,i,dt)+f2(s,i,dt)+f3(s,i,dt)+f4(s,i,dt));

%% Simulation
AllsimS=cell(NUMBER_OF_SIMUL,1); % Set the simulation 
AllsimI=cell(NUMBER_OF_SIMUL,1);
last_time_I=zeros(NUMBER_OF_SIMUL,1);
AllT=cell(NUMBER_OF_SIMUL,1);
for j = 1:NUMBER_OF_SIMUL
    clear T; clear simulS; clear simulI;
    simulS(1) = S0;
    simulI(1) = I0;
    i=2;
    T(1)=-log(rand)/(f1(S0,I0,1)+f2(S0,I0,1)+f3(S0,I0,1)+f4(S0,I0,1));
    while sum(T)<TOTAL_TIME
        pres=simulS(i-1); prei=simulI(i-1);
        c=rand;
        sumf = f1(pres,prei,T(i-1))+f2(pres,prei,T(i-1)) ...
            +f3(pres,prei,T(i-1))+f4(pres,prei,T(i-1));
        s1 = f1(pres,prei,T(i-1))/sumf;
        s2 = s1+f2(pres,prei,T(i-1))/sumf;
        s3 = s2+f3(pres,prei,T(i-1))/sumf;
        s4 = s3+f4(pres,prei,T(i-1))/sumf;
        if sumf == 0
            simulS(i) = pres;
            simulI(i) = prei;
        else
            if(c<=s1)
                simulS(i)=simulS(i-1)-1; simulI(i)=simulI(i-1)+1;
            elseif(c>s1 && c<=s2)
                simulS(i)=simulS(i-1); simulI(i)=simulI(i-1)-1;
            elseif(c>s2 && c<=s3)
                simulS(i)=simulS(i-1)+1; simulI(i)=simulI(i-1)-1;
            elseif(c>s3 && c<=s4)
                simulS(i)=simulS(i-1)+1; simulI(i)=simulI(i-1);
            end
        end
        
        sumf_next = f1(simulS(i),simulI(i),1)+f2(simulS(i),simulI(i),1)...
             +f3(simulS(i),simulI(i),1)+f4(simulS(i),simulI(i),1);
        if sumf_next == 0
            T(i)=0.01;
        else
            T(i)=-log(rand)/sumf_next;
        end
        i=i+1;
    end
    T=cumsum(T);
    T(length(T))=[];
    t=[0 T];
    AllsimS{j,1} = simulS;
    AllsimI{j,1} = simulI;
    last_time_I(j) = simulI(end);
    AllT{j,1} = t;
end

%% Quasi-stationary distribution (numerical)
distribution_of_I=zeros(N+1,1);
for i=0:N
    distribution_of_I(i+1)=sum(last_time_I==i);
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
figure(2);
hold on;                          % Figure 3.7
for j=1:NUMBER_OF_SIMUL
    plot(AllT{j,1},AllsimI{j,1});
end
plot(td,yd(:,2),'k--','LineWidth',2);
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
    cumDat = [cumDat AllsimI{i,1}];
end
h=histogram2(T, cumDat, 'BinWidth',[TIME_UNIT,INFECTIVE_UNIT],...
    'FaceColor','flat','DisplayStyle','bar3','ShowEmptyBins','on');
colorbar
plot(td,yd(:,2),'k-','LineWidth',4);
axis([0 X_AXIS 0 I_AXIS]);
xlabel('Time');
ylabel('Number of Infectives');
hold off;

% figure(4)
% hold on;
% t2=h.XBinEdges;
% t2(end)=[];
% Y=h.YBinEdges;
% Y(1)=[];
% count = h.Values;
% siz = size(count);
% maxData = max(count');
% for i= 1:siz(1)
%     count(i,:)=count(i,:)/maxData(i);
% end
% b=bar3(Y,count',1);
% for k=1:length(b)
%     zdata = b(k).ZData;
%     b(k).CData = zdata;
%     b(k).FaceColor = 'flat';
% end
% colorbar
% Xdat = get(b,'XData');
% for ii=1:length(Xdat)
%     Xdat{ii} = Xdat{ii}/(1/TIME_UNIT);
%     set(b(ii),'XData',Xdat{ii});
% end
% 
% plot(td,yd(:,2),'k-','LineWidth',4);
% axis([0 X_AXIS 0 I_AXIS]);
% xlabel('Time Step');
% ylabel('Number of Infectives');
% hold off;

%% Probability of an outbreak

I0_vals = [1,2,3];
kk=1;
prob_outbreak = zeros(3,1);
for I0 = I0_vals
    %% Simulation
    AllsimS=cell(NUMBER_OF_SIMUL,1); % Set the simulation 
    AllsimI=cell(NUMBER_OF_SIMUL,1);
    last_time_I=zeros(NUMBER_OF_SIMUL,1);
    AllT=cell(NUMBER_OF_SIMUL,1);
    for j = 1:NUMBER_OF_SIMUL
        clear T; clear simulS; clear simulI;
        simulS(1) = S0;
        simulI(1) = I0;
        i=2;
        T(1)=-log(rand)/(f1(S0,I0,1)+f2(S0,I0,1)+f3(S0,I0,1)+f4(S0,I0,1));
        while sum(T)<TOTAL_TIME
            pres=simulS(i-1); prei=simulI(i-1);
            c=rand;
            sumf = f1(pres,prei,T(i-1))+f2(pres,prei,T(i-1)) ...
                +f3(pres,prei,T(i-1))+f4(pres,prei,T(i-1));
            s1 = f1(pres,prei,T(i-1))/sumf;
            s2 = s1+f2(pres,prei,T(i-1))/sumf;
            s3 = s2+f3(pres,prei,T(i-1))/sumf;
            s4 = s3+f4(pres,prei,T(i-1))/sumf;
            if sumf == 0
                simulS(i) = pres;
                simulI(i) = prei;
            else
                if(c<=s1)
                    simulS(i)=simulS(i-1)-1; simulI(i)=simulI(i-1)+1;
                elseif(c>s1 && c<=s2)
                    simulS(i)=simulS(i-1); simulI(i)=simulI(i-1)-1;
                elseif(c>s2 && c<=s3)
                    simulS(i)=simulS(i-1)+1; simulI(i)=simulI(i-1)-1;
                elseif(c>s3 && c<=s4)
                    simulS(i)=simulS(i-1)+1; simulI(i)=simulI(i-1);
                end
            end

            sumf_next = f1(simulS(i),simulI(i),1)+f2(simulS(i),simulI(i),1)...
                 +f3(simulS(i),simulI(i),1)+f4(simulS(i),simulI(i),1);
            if sumf_next == 0
                T(i)=0.01;
            else
                T(i)=-log(rand)/sumf_next;
            end
            i=i+1;
        end
        T=cumsum(T);
        T(length(T))=[];
        t=[0 T];
        AllsimS{j,1} = simulS;
        AllsimI{j,1} = simulI;
        last_time_I(j) = simulI(end);
        AllT{j,1} = t;
    end

    %% Quasi-stationary distribution (numerical)
    distribution_of_I=zeros(N+1,1);
    for i=0:N
        distribution_of_I(i+1)=sum(last_time_I==i);
    end
    prob_outbreak(kk)=1-distribution_of_I(1)/NUMBER_OF_SIMUL;
    kk=kk+1;
end

figure(5);
bar([prob_outbreak, [1/2;(1-(1/2)^2);(1-(1/2)^3)]]);
xlabel("I_0"); ylabel("probability of an outbreak");
legend(["Simulation", "Theory"],"Location","northwest");
%%
I0=2;
R0=[0.5,2,5];
beta_vals = R0.*(bir+gam);
kk=1;
prob_outbreak=zeros(3,1);

for beta = beta_vals
    % Stochastic SIR
    f1 = @(s,i,dt) beta*i*s*dt/N;      % Set function handle for simplicity of coding
    f2 = @(s,i,dt) gam*i*dt;
    f3 = @(s,i,dt) bir*i*dt;
    f4 = @(s,i,dt) bir*(N-s-i)*dt;
    f5 = @(s,i,dt) 1-(f1(s,i,dt)+f2(s,i,dt)+f3(s,i,dt)+f4(s,i,dt));
    
    %% Simulation
    AllsimS=cell(NUMBER_OF_SIMUL,1); % Set the simulation 
    AllsimI=cell(NUMBER_OF_SIMUL,1);
    last_time_I=zeros(NUMBER_OF_SIMUL,1);
    AllT=cell(NUMBER_OF_SIMUL,1);
    for j = 1:NUMBER_OF_SIMUL
        clear T; clear simulS; clear simulI;
        simulS(1) = S0;
        simulI(1) = I0;
        i=2;
        T(1)=-log(rand)/(f1(S0,I0,1)+f2(S0,I0,1)+f3(S0,I0,1)+f4(S0,I0,1));
        while sum(T)<TOTAL_TIME
            pres=simulS(i-1); prei=simulI(i-1);
            c=rand;
            sumf = f1(pres,prei,T(i-1))+f2(pres,prei,T(i-1)) ...
                +f3(pres,prei,T(i-1))+f4(pres,prei,T(i-1));
            s1 = f1(pres,prei,T(i-1))/sumf;
            s2 = s1+f2(pres,prei,T(i-1))/sumf;
            s3 = s2+f3(pres,prei,T(i-1))/sumf;
            s4 = s3+f4(pres,prei,T(i-1))/sumf;
            if sumf == 0
                simulS(i) = pres;
                simulI(i) = prei;
            else
                if(c<=s1)
                    simulS(i)=simulS(i-1)-1; simulI(i)=simulI(i-1)+1;
                elseif(c>s1 && c<=s2)
                    simulS(i)=simulS(i-1); simulI(i)=simulI(i-1)-1;
                elseif(c>s2 && c<=s3)
                    simulS(i)=simulS(i-1)+1; simulI(i)=simulI(i-1)-1;
                elseif(c>s3 && c<=s4)
                    simulS(i)=simulS(i-1)+1; simulI(i)=simulI(i-1);
                end
            end

            sumf_next = f1(simulS(i),simulI(i),1)+f2(simulS(i),simulI(i),1)...
                 +f3(simulS(i),simulI(i),1)+f4(simulS(i),simulI(i),1);
            if sumf_next == 0
                T(i)=0.01;
            else
                T(i)=-log(rand)/sumf_next;
            end
            i=i+1;
        end
        T=cumsum(T);
        T(length(T))=[];
        t=[0 T];
        AllsimS{j,1} = simulS;
        AllsimI{j,1} = simulI;
        last_time_I(j) = simulI(end);
        AllT{j,1} = t;
    end

    %% Quasi-stationary distribution (numerical)
    distribution_of_I=zeros(N+1,1);
    for i=0:N
        distribution_of_I(i+1)=sum(last_time_I==i);
    end
    prob_outbreak(kk)=1-distribution_of_I(1)/NUMBER_OF_SIMUL;
    kk=kk+1;
end

figure(6);
bar([prob_outbreak, [0;(1-(1/2)^2);(1-(1/5)^3)]]);
xticklabels(R0);
xlabel("R_0"); ylabel("probability of an outbreak");
legend(["Simulation", "Theory"],"Location","northwest");