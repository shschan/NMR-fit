function [pfit, pfitErr, sim] = fit_peak_1state(w, y)

em = 5*pi; % line broadening factor 

% initial parameter estimates:
p0_1 = [ 0.9112   -1.2443   47.1003]; % amplitude, frequency (s-1), R2 (s-1) for first peak


% generic lorentzian:
L = @(A,w0,R2) (A*(R2+em).^2) ./ ((R2+em)^2 + (w-1000*w0).^2);

% fit three peaks
p0 = [p0_1 ];
ysim = @(p) L(p(1),p(2),p(3));  % -state

% p0 = [p0_1];
% ysim = @(p) L(p(1),p(2),p(3));   %1-state
resid = @(p) y - ysim(p);

% do the fitting
[pfit, ~, ~, ~, ~, ~, jac] = lsqnonlin(resid,p0,[],[],optimoptions(@lsqnonlin,'display','none'));

% calculate standard error
res = resid(pfit);
chi2 = sum(res.^2);
[rows,cols]=size(jac);
Sigma_o = sqrt((res'*res) / (rows-cols)); % (may need to switch multiplication order above depending if resid is row or column vector)
Q_xx = inv(jac'*jac);
pfitErr = full(Sigma_o .* sqrt(diag(Q_xx)));

sim = ysim(pfit);
