%%
clear all; clc; clf;

%% Set the Constants 
TOTAL_POPULATION=20; 
TRANSMISSIBILITY=1;   
BIRTH_RATE=0; 
RECOVERY_RATE=0.5; 
INTEREVENT_TIME=0.01;
TOTAL_TIME=100; 
NUMBER_OF_SIMUL=1000; 
ZERO_INFECTIVES=1;      
ZERO_SUSCEPTIBLES=TOTAL_POPULATION-1;

%% Drawing constants
X_AXIS = 20; I_AXIS = 50;

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

%% parameter changes
R0=[0.5,2,5];
beta_vals = R0.*(bir+gam);
kk=1;
final_sizes=zeros(N+1,3);

for beta = beta_vals
    %% Stochastic SIR
    f1 = @(s,i,dt) beta*i*s*dt/N;      % Set function handle for simplicity of coding
    f2 = @(s,i,dt) gam*i*dt;
    f3 = @(s,i,dt) bir*i*dt;
    f4 = @(s,i,dt) bir*(N-s-i)*dt;
    f5 = @(s,i,dt) 1-(f1(s,i,dt)+f2(s,i,dt)+f3(s,i,dt)+f4(s,i,dt));

    %% Simulation
    AllsimS=cell(NUMBER_OF_SIMUL,1); % Set the simulation 
    AllsimI=cell(NUMBER_OF_SIMUL,1);
    final_size=zeros(NUMBER_OF_SIMUL,1);
    AllT=cell(NUMBER_OF_SIMUL,1);
    for j = 1:NUMBER_OF_SIMUL
        clear T; clear simulS; clear simulI;
        simulS(1) = S0;
        simulI(1) = I0;
        i=2;
        num_of_cases = I0;
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
                    num_of_cases=num_of_cases+1;
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
        final_size(j) = num_of_cases;
        AllT{j,1} = t;
    end

    %% Final sizes (numerical)
    distribution_of_I=zeros(N+1,1);
    for i=0:N
        distribution_of_I(i+1)=sum(final_size==i);
    end
    final_sizes(:,kk)=distribution_of_I;
    kk=kk+1;
end
figure();
hold on;
plot(1:N, final_sizes(2:end,1)/NUMBER_OF_SIMUL, "ro--", "LineWidth", 2);
plot(1:N, final_sizes(2:end,2)/NUMBER_OF_SIMUL, "go--", "LineWidth", 2);
plot(1:N, final_sizes(2:end,3)/NUMBER_OF_SIMUL, "bo--", "LineWidth", 2);
xlabel("Final Size");
ylabel("Distribution");
legend(["R_0=0.5","R_0=2","R_0=5"]);