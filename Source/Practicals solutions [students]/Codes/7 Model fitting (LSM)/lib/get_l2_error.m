function squared_error = get_l2_error(seroprev_data, seroprev_est)
    squared_error = sqrt(sum((seroprev_data - seroprev_est).^2));
end