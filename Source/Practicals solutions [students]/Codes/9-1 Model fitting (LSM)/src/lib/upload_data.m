function output = upload_data(filename)
%upload_data - upload seroprevalence data
%
% Syntax: output = upload_data(filename)
%
% Input: filename(string)
% Output: data whose form is determined

%% Upload the data
data = csvread(filename);
size_check_pos = prod(size(data(:,1)) == size(data(:,2)));
size_check_tot = prod(size(data(:,1)) == size(data(:,3)));
assert(logical(size_check_pos),...
    "Data size error: data should have all the values");
assert(logical(size_check_tot),...
    "Data size error: data should have all the values");

%% Categorize the data
age = data(:,1);
num_of_positive = data(:,2);
num_of_samples = data(:,3);

%% Calculate the seroprevalence for each age
seroprev = num_of_positive ./ num_of_samples;

%% Output
output = [age, seroprev];
end