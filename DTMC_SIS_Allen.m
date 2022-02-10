% Discrete Time Markov Chain
% SIS Epidemic Model
% Transition Matrix and Graph of Probability Distribution
close all; clear all; clc
% set(gca,'FontSize',18);
% set(0,'DefaultAxesFontSize',18);
time=2000;
dtt=0.01; % Time step
beta=1*dtt;
b=0.25*dtt;
gama=0.25*dtt;
N=100; % Total population size
en=50; % plot every enth time interval
T=zeros(N+1,N+1); % T is the transition matrix, defined below
v=linspace(0,N,N+1);
p=zeros(time+1,N+1);
p(1,3)=1; % Two individuals initially infected.
bt=beta*v.*(N-v)/N;
dt=(b+gama)*v;
for i=2:N % Define the transition matrix
    T(i,i)=1-bt(i)-dt(i); % diagonal entries
    T(i,i+1)=dt(i+1); % superdiagonal entries
    T(i+1,i)=bt(i); % subdiagonal entries
end
T(1,1)=1;
T(1,2)=dt(2);
T(N+1,N+1)=1-dt(N+1);
for t=1:time
    y=T*p(t,:)';
    p(t+1,:)=y';
end
pm(1,:)=p(1,:);
for t=1:time/en
    pm(t+1,:)=p(en*t,:);
end
ti=linspace(0,time,time/en+1);
st=linspace(0,N,N+1);
mesh(st,ti,pm);
xlabel('Number of Infectives');
ylabel('Time Steps');
zlabel('Probability');
view(140,30);
axis([0,N,0,time,0,1]);