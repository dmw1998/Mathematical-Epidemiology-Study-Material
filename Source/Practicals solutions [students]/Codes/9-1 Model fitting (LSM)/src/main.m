%% make result path
mkdir("../res/");

%% clear all
clear; close all; clc;

%% add lib as path
addpath("./lib");

%% color definition (each plot)
clr.DATA_COLOR = [0.5, 0.5, 0.5];
clr.MODEL_COLOR = [0, 0.4470, 0.7410];
clr.INITIAL_COLOR = [0, 0.4470, 0.7410, 0.5];
clr.FIGURE_SIZE = [10, 10, 1360, 768];

%% Part 1
disp(" "); 
disp(" "); 
fprintf("Part 1 begins\n");

%% Q1
disp(" "); 

uk.data = upload_data(...
    "../data/Practical9-1 Model fitting-data/seroprevalence_uk.csv");
uk.ages = uk.data(:,1);
uk.sero_prev = uk.data(:,2);

uk.initial_guess = 0.12;
uk.initial_seroprev = seroprev(uk.ages, uk.initial_guess);

draw_seroprev(uk.ages, uk.sero_prev, uk.initial_seroprev);
fprintf("Q1. Squared error value is %g\n", ...
        sqval(uk.sero_prev, uk.initial_seroprev));
saveas(gcf, "../res/uk_initial_guess_seroprevalence.png", "png");
close gcf;

%% Q2
disp(" "); 

uk.est_foi = estimate_foi(uk.initial_guess, uk.data);   
fprintf("Q2. After estimation, FOI is %g\n", uk.est_foi);
uk.est_seroprev = seroprev(uk.ages, uk.est_foi);

draw_seroprev(uk.ages, uk.sero_prev, ...
        [uk.initial_seroprev, uk.est_seroprev]);
fprintf("Squared error value after calibration is %g\n", ...
        sqval(uk.sero_prev, uk.est_seroprev));
saveas(gcf, "../res/uk_best_fitting_seroprevalence.png", "png");
close gcf;

%% Q4
disp(" "); 

uk.A = 1/uk.est_foi;

fprintf("Q4. Average age at infection is %g\n", uk.A);

%% Q5
disp(" "); 

uk.L = 60;
uk.R0 = uk.L/uk.A;

fprintf("Q5. The basic reproduction number R_0 is %g\n", uk.R0);
fprintf("And the herd immunity threshold from it is %g %%\n", ...
    100*(1-1/uk.R0));

%% Q6
disp(" "); 

ch.data = upload_data(...
    "../data/Practical9-1 Model fitting-data/seroprevalence_china.csv");
ch.ages = ch.data(:,1);
ch.sero_prev = ch.data(:,2);

ch.initial_guess = 0.12;
ch.initial_seroprev = seroprev(ch.ages, ch.initial_guess);

fprintf("Q6. Squared error value is %g\n", ...
        sqval(ch.sero_prev, ch.initial_seroprev));
    
ch.est_foi = estimate_foi(ch.initial_guess, ch.data); 

fprintf("After estimation, FOI is %g\n", ch.est_foi);

ch.est_seroprev = seroprev(ch.ages, ch.est_foi);

draw_seroprev(ch.ages, ch.sero_prev, ...
        [ch.initial_seroprev, ch.est_seroprev]);
fprintf("Squared error value after calibration is %g\n", ...
        sqval(ch.sero_prev, ch.est_seroprev));
saveas(gcf, "../res/ch_best_fitting_seroprevalence.png", "png");
close gcf;

ch.A = 1/ch.est_foi;

fprintf("Average age at infection is %g\n", ch.A);

ch.L = uk.L; % Assumption
ch.R0 = ch.L/ch.A;

fprintf("The basic reproduction number R_0 is %g\n", ch.R0);
fprintf("And the herd immunity threshold from it is %g %%\n", ...
    100*(1-1/ch.R0));

%% Q8,9
disp(" "); 

uk.initial_seroprev_ma = seroprev_ma(uk.ages, uk.initial_guess);
uk.est_foi_ma = estimate_foi_ma(uk.initial_guess, uk.data);

fprintf("Q8,9[UK]. After estimation with maternal immunity, FOI is %g\n", ...
        uk.est_foi_ma);
    
uk.est_seroprev_ma = seroprev_ma(uk.ages, uk.est_foi);

draw_seroprev(uk.ages, uk.sero_prev, ...
        [uk.initial_seroprev, uk.est_seroprev_ma]);
fprintf("Squared error value after calibration is %g\n", ...
        sqval(uk.sero_prev, uk.est_seroprev_ma));
saveas(gcf, "../res/uk_best_fitting_seroprevalence_ma.png", "png");
close gcf;

ch.initial_seroprev_ma = seroprev_ma(ch.ages, ch.initial_guess);
ch.est_foi_ma = estimate_foi_ma(ch.initial_guess, ch.data);   

fprintf("Q8,9[China]. After estimation with maternal immunity, FOI is %g\n", ...
        ch.est_foi_ma);
    
ch.est_seroprev_ma = seroprev_ma(ch.ages, ch.est_foi);

draw_seroprev(ch.ages, ch.sero_prev, ...
        [ch.initial_seroprev, ch.est_seroprev_ma]);
fprintf("Squared error value after calibration is %g\n", ...
        sqval(ch.sero_prev, ch.est_seroprev_ma));
saveas(gcf, "../res/ch_best_fitting_seroprevalence_ma.png", "png");
close gcf;

%% Q10
disp(" "); 
disp("Q10.");

uk.neg_log_suscep_tot = -log(1-uk.sero_prev); % Susceptible proportion
ch.neg_log_suscep_tot = -log(1-ch.sero_prev);

fig = figure("Units","pixels","Position",clr.FIGURE_SIZE);
plt_data = plot(uk.ages, uk.neg_log_suscep_tot, ...
    "o", "MarkerSize", 5, "LineWidth", 2, "Color", clr.DATA_COLOR);
set(plt_data, "markerfacecolor", get(plt_data, "color")); 

xlabel("Ages");
ylabel("-log(S_a/N_a)");
set(findall(fig,"-property","FontSize"),"FontSize",20);

saveas(fig, "../res/uk_neg_log_suscept_tot.png", "png");
clear fig;
close all;

fig = figure("Units","pixels","Position",clr.FIGURE_SIZE);
plt_data = plot(ch.ages, ch.neg_log_suscep_tot, ...
    "o", "MarkerSize", 5, "LineWidth", 2, "Color", [0.5, 0.5, 0.5]);
set(plt_data, "markerfacecolor", get(plt_data, "color")); 

xlabel("Ages");
ylabel("-log(S_a/N_a)");
set(findall(fig,"-property","FontSize"),"FontSize",20);

saveas(fig, "../res/ch_neg_log_suscept_tot.png", "png");
clear fig;
close all;

%% Q11
disp(" "); 

uk.initial_guess = [0.12, 0.12];
uk.initial_seroprev_diff = seroprev(uk.ages, uk.initial_guess, 15);

fprintf("Q11. Squared error value before calibration is %g\n", ...
        sqval(uk.sero_prev, uk.initial_seroprev_diff));
    
uk.est_foi_diff = estimate_foi(uk.initial_guess, uk.data, 15); 

fprintf("After estimation, FOI before 15 is %g\n", uk.est_foi_diff(1));
fprintf("FOI after 15 is %g\n", uk.est_foi_diff(2));
uk.est_seroprev_diff = seroprev(uk.ages, uk.est_foi_diff,15);

draw_seroprev(uk.ages, uk.sero_prev, ...
        [uk.initial_seroprev_diff, uk.est_seroprev_diff]);
fprintf("Squared error value after calibration is %g\n", ...
        sqval(uk.sero_prev, uk.est_seroprev_diff));
saveas(gcf, "../res/uk_best_fitting_seroprevalence_2.png", "png");

disp(" ");
ch.initial_guess = [0.12, 0.12];
ch.initial_seroprev_diff = seroprev(ch.ages, ch.initial_guess,15);

fprintf("Squared error value before calibration is %g\n", ...
        sqval(ch.sero_prev, ch.initial_seroprev_diff));

ch.est_foi_diff = estimate_foi(ch.initial_guess, ch.data, 15);   

fprintf("After estimation, FOI before 15 is %g\n", ch.est_foi_diff(1));
fprintf("FOI after 15 is %g\n", ch.est_foi_diff(2));

ch.est_seroprev_diff = seroprev(ch.ages, ch.est_foi_diff,15);

draw_seroprev(ch.ages, ch.sero_prev, ...
        [ch.initial_seroprev_diff, ch.est_seroprev_diff]);
fprintf("Squared error value after calibration is %g\n", ...
        sqval(ch.sero_prev, ch.est_seroprev_diff));
saveas(gcf, "../res/ch_best_fitting_seroprevalence_2.png", "png");
close all;

%%

%% Part 2
disp(" "); 
disp(" "); 

fprintf("Part 2 begins\n");

p2.data = upload_data2( ...
    "../data/Practical9-1 Model fitting-data/incidence_measles.csv");
p2.params = set_params();
p2.initial_guess = 1e-4;
p2.params.beta = optimize_beta(p2.initial_guess, p2.params, p2.data);

fprintf("Q1. The best-fitting value of transmission rate is %g\n",...
    p2.params.beta);

p2.init_sol = solve_SEIR(p2.initial_guess,p2.params,p2.data(end,1));
p2.init_inc = simul_inc(p2.init_sol,p2.params.f);

p2.opt_sol = solve_SEIR(p2.params.beta,p2.params,p2.data(end,1));
p2.opt_inc = simul_inc(p2.opt_sol,p2.params.f);

%% drawing

figure("Units", "pixels", "Position", clr.FIGURE_SIZE);
hold on;
plt_data = plot(p2.data(:,1), p2.data(:,2), ...
    "o", "MarkerSize", 5, "LineWidth", 2, "Color", clr.DATA_COLOR);
set(plt_data, "markerfacecolor", get(plt_data, "color")); 
plot(p2.data(:,1), p2.init_inc, ...
    "LineWidth", 2, "Color", clr.INITIAL_COLOR);
plot(p2.data(:,1), p2.opt_inc, ...
    "LineWidth", 2, "Color", clr.MODEL_COLOR);

xlabel("Time [Days]");
ylabel("Weekly incidence of measles");
legend(["Data", "Init.", "Est."], "Location", "northeastoutside");
set(findall(gcf,"-property","FontSize"),"FontSize",20);

saveas(gcf, "../res/fitting_incidence_with_LSM.png", "png");
close all;
clear clr;
clear plt_data;