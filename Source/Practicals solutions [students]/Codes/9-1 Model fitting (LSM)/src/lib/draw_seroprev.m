function draw_seroprev(ages, seroprev_data, seroprev_est)
%draw_seroprev - draw seroprevalence by age from data and estimated FOI 
%
% Syntax: draw_seroprev(ages, seroprev_data, seroprev_est)
%
% Input: ages, seroprev_data, seroprev_est (all should have the same size.)
DATA_COLOR = [0.5, 0.5, 0.5];
MODEL_COLOR = [0, 0.4470, 0.7410];
INITIAL_COLOR = [0, 0.4470, 0.7410, 0.5];

figure("Units","pixels","Position",[10,10,1360,768]);
hold on;

plt_data = plot(ages, 100*seroprev_data, "o", ...
        "MarkerSize", 5, "LineWidth", 2, "Color", DATA_COLOR);
set(plt_data, "markerfacecolor", get(plt_data, "color")); 

switch size(seroprev_est, 2)
    case 1
        plot(ages, 100*seroprev_est, "-", ...
            "LineWidth", 2, "Color", MODEL_COLOR);
        legend_str = ["Data", "Model"];
    case 2
        plot(ages, 100*seroprev_est(:,1), "-", ...
            "LineWidth", 2, "Color", INITIAL_COLOR);
        plot(ages, 100*seroprev_est(:,2), "-", ...
            "LineWidth", 2, "Color", MODEL_COLOR);
        legend_str = ["Data", "Before Est.", "Est."];
    otherwise
        error("Drawing Error: For drawing seroprevalence from estimated FOI, the input should have less than or equal to 2 column number.")
end

xlim([-0.5, 45])
ylim([-0.5,100.5]);
xlabel("Ages");
ylabel("Seroprevalence [%]");
legend(legend_str, "Location", "northeastoutside");
set(findall(gcf,"-property","FontSize"),"FontSize",20);
end