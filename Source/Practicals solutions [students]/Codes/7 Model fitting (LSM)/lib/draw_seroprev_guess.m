function draw_seroprev_guess(ages, seroprev_data, seroprev_guess)
    % set the plot color, never mind if you don't want to use these.
    DATA_COLOR = [0.5, 0.5, 0.5]; % grey
    INITIAL_COLOR = [0, 0.4470, 0.7410]; % blue
    
    figure("Units","pixels","Position",[10,10,640,360]);
    hold on;
    
    plot(ages, 100*seroprev_data, "o", ...
        "MarkerSize", 5, "LineWidth", 2, "Color", DATA_COLOR);
    
    plot(ages, 100*seroprev_guess, "-", ...
        "LineWidth", 2, "Color", INITIAL_COLOR);
    legend_str = ["Data", "Initial Guess"];
    
    xlim([-0.5, 45])
    ylim([-0.5, 101]);
    xlabel("Ages");
    ylabel("Seroprevalence [%]");
    legend(legend_str, "Location", "northeastoutside");
    set(findall(gcf,"-property","FontSize"),"FontSize",8);
end