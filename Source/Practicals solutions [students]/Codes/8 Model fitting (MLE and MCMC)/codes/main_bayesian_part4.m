clc; clear all; close all;

%% PART 2 - measles incidence

% data
T = readmatrix('incidence_measles.csv');
time_stamp = T(:,1);        % data time
time_stampc = 0:140;   % continuous time stamp for prediction
data = T(:,2);              % positive

param.time_stamp = time_stamp;
param.time_stampc = time_stampc;
param.data = data;

% fixed parameter
N = 1e+5;
initial = [N-1 0 1 0];
f = 1/8;
r = 1/7;

param.initial = initial;
param.f = f;
param.r = r;

%% Question 13-1
% 13-1-1 MCMC procedure
f_mh = @(theta) sum(log(poisspdf(data,solve_SEIR(theta,param))));
N = 1e+4;   % number of samples
ini = 3e-5*rand(1,10);   % initails

for k = 1 : 10
    burn_in = 0;            % no burn-in
    NN = N + burn_in;
    C0 = 1e-15;
    candidate.mean = zeros(1,1);
    candidate.cov = C0;
    THETA = ini(k);
    rng('default')
    fprintf('start MCMC procedure \n')
    tic;
    for i = 2 : NN
        % draw a theta_prime from candidate distribution
        theta_prime = -ones(1,1);
        while theta_prime(1)<0
            theta_prime = mvnrnd(candidate.mean + THETA(:,i-1) , candidate.cov)';
        end
        % calculate Alpha
        % step 1 : target
        temp_target = exp(f_mh(theta_prime)-f_mh(THETA(:,i-1)));
        % step 2 : candidate
        temp_candidate = 1;     % random-walk
        
        % acceptance probability
        A = min ( 1, temp_target * temp_candidate );
        
        % draw a sample from uniform dist.
        uni = rand(1);
        
        % accept or reject
        if uni < A
            THETA(:,i) = theta_prime;
        else
            THETA(:,i) = THETA(:,i-1);
        end
    end
    t_mh = toc;
    fprintf('Computing time: %f seconds\n', t_mh)
    samples(k,:) = THETA(:,burn_in+1:end);
end

    % 13-1-2 Trace and Autocorrelation
    figure1 = figure('pos',[10 10 600 400]);
    plot(1:length(samples),samples)
    grid on; grid minor;
    ylabel('Samples')
    title('Trace')
    set(gca, 'FontSize', 15)
    saveas(gca, 'Q13-1 trace_for_various_initials','epsc')


%% Question 13-2
% 13-2-1 MCMC procedure
fprintf('start MLE procedure \n')
f_mh = @(theta) sum(log(poisspdf(data,solve_SEIR(theta,param))));
N = [1e+2 1e+3 1e+4 1e+5];  % number of samples

for k = 1 : 4
    burn_in = N(k)*0.05;        % 5% burn-in
    NN = N(k) + burn_in;
    C0 = 1e-15;
    candidate.mean = zeros(1,1);
    candidate.cov = C0;
    THETA = 1.8e-5;
    rng('default')
    fprintf('start MCMC procedure \n')
    tic;
    for i = 2 : NN
        % draw a theta_prime from candidate distribution
        theta_prime = -ones(1,1);
        while theta_prime(1)<0
            theta_prime = mvnrnd(candidate.mean + THETA(:,i-1) , candidate.cov)';
        end
        % calculate Alpha
        % step 1 : target
        temp_target = exp(f_mh(theta_prime)-f_mh(THETA(:,i-1)));
        % step 2 : candidate
        temp_candidate = 1;     % random-walk
        
        % acceptance probability
        A = min ( 1, temp_target * temp_candidate );
        
        % draw a sample from uniform dist.
        uni = rand(1);
        
        % accept or reject
        if uni < A
            THETA(:,i) = theta_prime;
        else
            THETA(:,i) = THETA(:,i-1);
        end
    end
    t_mh = toc;
    fprintf('Computing time: %f seconds\n', t_mh)
    
    samples = THETA(:,burn_in+1:end);
       
    % 13-2 Trace and Autocorrelation
    figure1 = figure('pos',[10 10 1200 400]);
    subplot(1,2,1)
    plot(1:length(samples),samples)
    grid on; grid minor;
    ylabel('Samples')
    title('Trace')
    set(gca, 'FontSize', 15)
    subplot(1,2,2)
    autocorr(samples, N(k)*0.05)
    set(gca, 'FontSize', 15)
    filename = sprintf('%3.0e',N(k));
    saveas(gca, sprintf('Q13-2 trace_autocorre_%s',filename),'epsc')
    
    % 13-3 distribution
    t = sort(samples);
    for i = 1 : length(t)
        yt(i) = f_mh(t(i));
    end
    figure1 = figure('pos',[10 10 600 400]);
    [counts,centers] = hist(samples,20);
    bar(centers,counts,0.8); hold on;
    plot(t,exp(yt-max(yt))*max(counts),'LineWidth',2);
    grid on; grid minor;
    ylabel('Samples')
    title('Histogram and Posterior distribution')
    set(gca, 'FontSize', 15)
    filename = sprintf('%3.0e',N(k));
    saveas(gca, sprintf('Q13-3 hist_distribution_%s',filename),'epsc')
end


%% Question 14

% 14-1 best fitting values - foi, average age of infection, R0, H
theta_MCMC = mean(samples');    % sample mean
R0 = theta_MCMC*initial(1)/r;    % basic reproduction number
H = 1-1/R0;                     % herd immunity thershold

fprintf('Q14. best-fitting value : %e \n',theta_MCMC);
fprintf('Q14. basic reproduction number : %f \n',R0)
fprintf('Q14. herd immunity threshold   : %f \n\n',H)

% 14-2 95% Credible interval
temp = sort(samples');
theta_MCMC_ci = [temp(ceil(length(temp')*0.025),:)' temp(floor(length(temp')*0.975),:)'];
R0_ci = theta_MCMC_ci*initial(1)./r;
H_ci = 1-1./R0_ci;

fprintf('Q14. CI for force of infection : [%e %e] \n',theta_MCMC_ci(1),theta_MCMC_ci(2));
fprintf('Q14. CI for basic reproduction number : [%f %f] \n',R0_ci(1),R0_ci(2))
fprintf('Q14. CI for herd immunity threshold   : [%f %f] \n\n',H_ci(1),H_ci(2))

% 14-3 best-fitting plot
y_theta_MCMC = solve_SEIR(theta_MCMC,param);

figure1 = figure('pos',[10 10 600 400]);
plot(time_stamp,data,'.','MarkerSize',20);  % data plot
hold on; 
plot(time_stamp,y_theta_MCMC,'LineWidth',2)    % prediction plot
temp_legend = sprintf('Prediction (\\beta = %e)',theta_MCMC);
legend('Data',temp_legend)
xlabel('Time (day)')
ylabel('Weekly incidence of measles')
grid on; grid minor;
set(gca, 'FontSize', 15)
saveas(gca, 'Q14-1 fitting', 'epsc')

