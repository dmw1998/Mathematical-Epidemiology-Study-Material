function squared_error = get_l2_error_from_foi(ages, ...
                                                    seroprev_data, ...
                                                    foi)
    seroprev_from_foi = get_seroprev_from_foi(ages, foi);
    squared_error = get_l2_error(seroprev_from_foi, seroprev_data);
end