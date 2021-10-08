function seroprev = get_seroprev_from_foi_with_maternality(age, foi) 
    seroprev = 1 - exp(- foi .* max(age-0.5, 0));
end