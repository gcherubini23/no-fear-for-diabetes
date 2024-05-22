function tau = delay(s1, s2)
    [r, lags] = xcorr(s1 - mean(s1), s2 - mean(s2), 'coeff');
    [~, I] = max(abs(r));
    % figure
    % plot(lags, r, '-')
    tau = lags(I);
end

