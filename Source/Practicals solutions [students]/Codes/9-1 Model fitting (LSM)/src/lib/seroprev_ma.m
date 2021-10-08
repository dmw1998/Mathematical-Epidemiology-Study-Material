function output = seroprev_ma(ages, foi, diff_age)
%seroprev - calculate seroprevalence for each age with given FOI assuming
%           maternal immunity.
%
% Syntax: output = seroprev_ma(ages, foi, diff_age)
%
% Input: given ages and foi (double), diff_age means the age values whose FOI is differentiated when foi is not a scalar.
% Output: seroprevalence value for each age with given foi. [sero]

%% Check error.
if nargin > 2
    check_size_diff_age = prod(size(foi)-[0,1] == size(diff_age));
    assert(logical(check_size_diff_age), ...
        "Given FOI Error: The ages differentiating foi values should have the one less size of the given foi.")
else
    diff_age = ages(end);
end

%% Calcualte the foi
sero = zeros(size(ages));

check_points = find(ages-diff_age>0,1);
if isempty(check_points)
    check_points = length(ages)+1;
end
sero(1:check_points-1) = 1-exp(-foi(1)*(ages(1:check_points-1)-0.5));

if (check_points~=length(ages)+1)
    sero(check_points:end) = ...
        1-exp(-diff_age*(foi(1)-foi(2)))*exp(-foi(2)*(ages(check_points:end)-0.5));
end

%% output the data
output = sero;
end