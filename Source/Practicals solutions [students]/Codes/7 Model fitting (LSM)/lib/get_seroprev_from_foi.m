function seroprev = get_seroprev_from_foi(age, foi) 
    seroprev = 1 - exp(- foi .* age);
end