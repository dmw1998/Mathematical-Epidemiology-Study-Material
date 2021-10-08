function estimated_foi = estimate_foi_ma(initial_guess, data, diff_age)
%estimate_foi - estimate the foi values from data initial guess assuming
%               the perfect maternal immunity
%
% Syntax: estimated_foi = estimate_foi_ma(initial_guess, data)
%
% Input: initial_guess (double), data(vector)
% Output: Estimated FOI value who minimizes negative log-likelihood of data.

ages = data(:,1);
sero_prev = data(:,2);
if nargin == 2
    estimated_foi = fminsearchbnd(@(x) calculate_sq(x, ages, sero_prev), ...
                                initial_guess, ...
                                zeros(size(initial_guess)) ...
                                );
else
    estimated_foi = fminsearchbnd(@(x) calculate_sq(x, ages, sero_prev, diff_age), ...
                                initial_guess, ...
                                zeros(size(initial_guess)) ...
                                );
end
end

function sq = calculate_sq(guessed_foi, ages, sero_prev, diff_age)
if nargin == 3
    output = seroprev_ma(ages, guessed_foi);
else
    output = seroprev_ma(ages, guessed_foi, diff_age);
end
sq = sqval(sero_prev, output);
end