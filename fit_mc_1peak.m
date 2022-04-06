function [pfit, pfitErr, sim0] = fit_mc_1peak(w, y)
% w = frequency (in s-1), y = signal
% re-center w to be near zero
% scale down signal intensities (y) so they're around 1

N_repeats = 200;

do_plot = true;

if do_plot
    figure
    plot(w,0*y,'-','color',[.7 .7 .7])
    hold on
    plot(w, y, 'kx')
end

% do first fit to raw data
[pfit0, pfitErr0, sim0] = fit_peak_1state(w, y);

resid0 = y - sim0;

N = length(w);
block_size = 400;
N_blocks = N - block_size;

for j=1:N_repeats
    new_y = sim0;
    for i=1:ceil(N/block_size)
        idx = 1+(i-1)*block_size:min(N,i*block_size);
        block_idx = randi(N_blocks);
        dy = resid0(block_idx:block_idx+block_size-1);
        dy = dy(1:length(idx));
        new_y(idx) = new_y(idx) + dy;
    end
    %plot(w, new_y, 'cx')
    [pfitMC{j}, ~, sim] = fit_peak_1state(w, new_y);
    if do_plot
        %plot(w, new_y, 'r-')
        plot(w, sim, 'b-')
    end
end
if do_plot
    plot(w, sim0, 'r-')
    hold off
    set(gca,'xdir','reverse');
    xlabel('frequency / s^{-1}')
    xlim([w(end) w(1)])
end

pfit = mean(cell2mat(pfitMC'));
pfitErr = std(cell2mat(pfitMC'));
