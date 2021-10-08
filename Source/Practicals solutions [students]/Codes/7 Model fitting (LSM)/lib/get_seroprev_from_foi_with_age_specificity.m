function seroprev = get_seroprev_from_foi_with_age_specificity(age, foi) 
    ex = - [foi(1)*age(age<15); 
            15*(foi(1)-foi(2))+foi(2)*age(age>=15)];
    seroprev = 1 - exp(ex);
end