function output = upload_data2(filename)
%upload_data - upload incidence data
%
% Syntax: output = upload_data2(filename)
%
% Input: filename(string)
% Output: data whose form is determined

%% Upload the data
data = csvread(filename);
size_check_inc = prod(size(data(:,1)) == size(data(:,2)));
assert(logical(size_check_inc),...
    "Data size error: data should have all the values");

%% Categorize the data
age = data(:,1);
incidence = data(:,2);

%% Output
output = [age, incidence];
end