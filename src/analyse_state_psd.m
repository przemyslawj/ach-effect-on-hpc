for i = 1:nexp
    if mod(i,3) < 2
        continue
    end
    
    % light off at experiment i - 1, light on at i
    pow_theta_light_off = psd_table.pow_theta(i-1,:);
    pow_theta_light_on = psd_table.pow_theta(i,:);
    [h, pval] = ttest(pow_theta_light_off, pow_theta_light_on);
end

groups = [ zeros(1, numel(pow_theta_light_off)) ones(1, numel(pow_theta_light_off))];
figure;
gscatter(groups, [pow_theta_light_off pow_theta_light_on], groups);