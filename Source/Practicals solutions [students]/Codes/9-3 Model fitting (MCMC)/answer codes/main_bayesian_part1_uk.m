clc; clear all; close all;

%% PART 1 - Seroprevalence UK

% data
T = readmatrix('seroprevalence_uk.csv');
time_stamp = T(:,1);        % data age
time_stampc = 0:0.01:44;    % continuous time stamp for prediction
Pa = T(:,2);                % positive
Na = T(:,3);                % N
data = Pa./Na;              % seroprevalence

%% Question 1
% 1-1 MCMC procedure
fprintf('start MLE procedure \n')
f_mh = @(theta) sum(log(binopdf(Pa,Na,1-exp(-theta.*time_stamp))));
N = [1e+2 1e+3 1e+4 1e+5];  % number of samples

for k = 1 : 1
    burn_in = N(k)*0.05;        % 5% burn-in
    NN = N(k) + burn_in;
    C0 = 8e-6;
    candidate.mean = zeros(1,1);
    candidate.cov = C0;
    THETA = 0.1;
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
    t_mh = toc
    fprintf('Computing time: %f seconds\n', t_mh)
    
    samples = THETA(:,burn_in+1:end);
       
    % 1-2 Trace and Autocorrelation
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
    saveas(gca, sprintf('Q1-1 trace_autocorre_%s',filename),'epsc')
    
    % 1-3 distribution
    t = sort(samples);
    for i = 1 : length(t)
        yt(i) = exp(f_mh(t(i)));
    end
    figure1 = figure('pos',[10 10 600 400]);
    [counts,centers] = hist(samples,20);
    bar(centers,counts,0.8); hold on;
    plot(t,yt/max(yt)*max(counts),'LineWidth',2);
    grid on; grid minor;
    ylabel('Samples')
    title('Histogram and Posterior distribution')
    set(gca, 'FontSize', 15)
    filename = sprintf('%3.0e',N(k));
    saveas(gca, sprintf('Q1-2 hist_distribution_%s',filename),'epsc')
    
end
%% Question 2

% 2-1 best fitting values - foi, average age of infection, R0, H
theta_MCMC = mean(samples');    % sample mean
A = 1/theta_MCMC;                % average age of infection
L = 60;                         % life expentancy
R0 = L/A;                       % basic reproduction number
H = 1-1/R0;                     % herd immunity thershold

fprintf('Q2. best-fitting value : %f \n',theta_MCMC);
fprintf('Q2. average age of infection  : %f \n',A)
fprintf('Q2. basic reproduction number : %f \n',R0)
fprintf('Q2. herd immunity threshold   : %f \n\n',H)

% 2-2 95% Credible interval
temp = sort(samples');
theta_MCMC_ci = [temp(ceil(length(temp')*0.025),:)' temp(floor(length(temp')*0.975),:)'];
A_ci = 1./theta_MCMC_ci;
R0_ci = L./A_ci;
H_ci = 1-1./R0_ci;

fprintf('Q2. CI for force of infection : [%f %f] \n',theta_MCMC_ci(1),theta_MCMC_ci(2));
fprintf('Q2. CI for average age of infection  : [%f %f] \n',A_ci(2),A_ci(1))
fprintf('Q2. CI for basic reproduction number : [%f %f] \n',R0_ci(1),R0_ci(2))
fprintf('Q2. CI for herd immunity threshold   : [%f %f] \n\n',H_ci(1),H_ci(2))

% 2-3 best-fitting plot
y_theta_MCMC = solve_catalytic(theta_MCMC,time_stampc);

figure1 = figure('pos',[10 10 600 400]);
plot(time_stamp,data,'.','MarkerSize',20);  % data plot
hold on; 
plot(time_stampc,y_theta_MCMC,'LineWidth',2)    % prediction plot
temp_legend = sprintf('Prediction (\\lambda = %f)',theta_MCMC);
legend('Data',temp_legend,'Location','east')
xlabel('Age (yrs)')
ylabel('Proportion positive')
grid on; grid minor;
set(gca, 'FontSize', 15)
saveas(gca, 'Q2-1 fitting', 'epsc')






