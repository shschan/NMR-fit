
% To baseline and fit a 1D 19F pulse-acquire NMR spectrum
% For more details, see Chan et al, Nat Chem 2022

load data.txt;   % process test.fid using nmrproc.com (ensure baselining at e.g. -70 and -50 ppm) and 
                %  convert test.ft1 file to data.txt using the following command: pipe2txt.tcl -index ppm test.ft1>data.txt


NS = 10* 1024;   % number of experiments summed times the NS for each experiment
conc = 10;       % conc in uM

normalise = NS * conc;

%%

% Quick look at spectrum, and select regions to fit baseline correction
figure
plot(data(:,1),data(:,2));
hold on

baseline_raw = data;                            %select region with signal
baseline_raw(5500:9600,:)=[];
plot(baseline_raw(:,1),baseline_raw(:,2),'r.'); %check signal is not selected

%%

% Baseline using 4th order polynomial

% Fit: 'untitled fit 1'.

xData = baseline_raw(:,1);
yData = baseline_raw(:,2);

syms x y

% Set up fittype and options.e

ft = fittype( 'poly4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

[fitresult, gof] = fit( xData, yData, ft, opts );


fitvalues = coeffvalues(fitresult);         %results of fit

baseline_fit = data;
data_corr = data;
for i=1:length(baseline_fit(:,1)),
    baseline_fit(i,2) = fitvalues(1)*baseline_fit(i,1)^4+ fitvalues(2)*baseline_fit(i,1)^3 + fitvalues(3)*baseline_fit(i,1)^2 + fitvalues(4)*baseline_fit(i,1) +fitvalues(5) ;
    data_corr(i,2) = data(i,2) - baseline_fit(i,2);
end

figure
subplot(3,1,1)
plot(data(:,1),data(:,2),'k-')
set(gca,'xdir','reverse');
    xlim([-70 -50])
    set(gca,'ytick',[])
    set(gca,'tickdir','out')
    set(gca,'box','off')
    hold off
subplot(3,1,2)
% plot(baseline_raw(:,1),baseline_raw(:,2),'k-')
plot(data(:,1),data(:,2),'k-')
hold on
plot(baseline_fit(:,1),baseline_fit(:,2),'r-')
set(gca,'xdir','reverse');
    xlim([-70 -50])
    set(gca,'ytick',[])
    set(gca,'tickdir','out')
    set(gca,'box','off')
    hold off
subplot(3,1,3)
hold on
plot(data_corr(:,1),data_corr(:,2),'k-')
set(gca,'xdir','reverse');
syms x
y= x*0;
plot(data_corr(:,1),data_corr(:,2),'k-')
fplot(y)
set(gca,'xdir','reverse');
    xlim([-70 -50])
    set(gca,'ytick',[])
    set(gca,'tickdir','out')
    set(gca,'box','off')
    hold off

data_extract = [data_corr(:,1) data_corr(:,2)./normalise];

% Adjust phasing of spectrum and repeat if needed

%%

% Fitting of spectrum for single state

data=data_extract;
clear xData
clear yData
    
for i=1:length(data_extract(:,1)),
    xData(i,1) = data(i,1)+62.2; %shift to set centre between 2 peaks as 0 Hz (change 62.2 to other if needed (and below too then))
 
end

yData = data(:,2)./max(data(:,2));

syms x y
x1 = -4;
x2 = +4;

% Set up fittype and options. 
ft = fittype( 'h1/(1+((CS1-x)/LW1)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

opts.StartPoint = [-0.4 0.1 1]; %adjust starting parameters here [CS(this is CS+62.2 as above) LW height]

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
fitvalues = coeffvalues(fitresult);         %results of fit

% Plot fit with data.
figure
subplot(2,1,1)
hold on
h = plot(xData-62.2, yData.*max(data(:,2)));
f1=fitvalues(3).*max(data(:,2))/(1+((fitvalues(1)-62.2-x)/fitvalues(2))^2);
fplot(f1,'b','linewidth',0.5)
set(gca,'xdir','reverse');
xlim([-64 -60])
grid on


%residuals plot
for i=1:length(yData),
    fitted(i,1) = fitvalues(3)/(1+((fitvalues(1)-xData(i))/fitvalues(2))^2); 
    residual(i,1) = yData(i,1)-fitted(i,1);
end
subplot(2,1,2)
hold on
h = plot(xData-62.2, residual.*max(data(:,2)),'color',[0.1 0.1 0.1].*5);
title('residuals')
set(gca,'xdir','reverse');
xlim([-64 -60])
    set(gca,'tickdir','out')
    set(gca,'box','off')
% Label axes
xlabel x
ylabel y
grid on


clear fiterror

confidence = confint(fitresult);
[m n] = size(confidence);
for i=1:n,
    fiterror(1,i) = abs(confidence(2,i)-confidence(1,i))/2;       %error of fit values
end

fit_data = fitvalues';
fit_data_err =  fiterror';

CS = [fit_data(1)-62.2 fit_data_err(1)]
LW = abs([fit_data(2) fit_data_err(2)].*470.611*2)
h = [fit_data(3) fit_data_err(3)].*max(data(:,2));

for i=1:1,
    integral(i,1) = LW(i)*h(i);
end

results = [CS(:,1) LW(:,1) populations_frac]

%%



%error analysis by bootstrapping

%initial parameters for fit - to be hardcoded into fit_peak_1state.m
p0_1 = [fit_data(3) fit_data(1).*470.611*2*pi/1000 abs(fit_data(2)*470.611*2*pi)]

xData_bs = xData.*470.611*2*pi; %convert ppm to s-1
yData_bs = yData;  

clear pfit pfitErr
[pfit, pfitErr] = fit_mc_1peak(xData_bs, yData_bs);

%%

fit_bs = reshape(pfit, [3 1])';
fit_bs_err = reshape(pfitErr, [3 1])';

%convert to units
h_bs = [fit_bs(:,1).*max(data(:,2)) fit_bs_err(:,1).*max(data(:,2))]
CS_bs = [fit_bs(:,2)*1000./(470.611*2*pi)-62.2 fit_bs_err(:,2)*1000./(470.611*2*pi)]
LW_bs = [fit_bs(:,3)./pi+10 fit_bs_err(:,3)./pi]

for i=1,
    integral_bs(i,1) = LW_bs(i,1)*h(i,1);
    integral_bs_err(i,1) = ((LW_bs(i,2)/LW_bs(i,1))^2 + (h_bs(i,2)/h_bs(i,1))^2)^0.5 * integral_bs(i,1);
    integral_sqerror(i,1) = (integral_bs_err(i)/integral_bs(i))^2;
end

populations_frac_bs = integral_bs./sum(integral_bs);

for i=1
    populations_frac_bs_err(i,1) = populations_frac_bs(i)*((integral_bs_err(i)/integral_bs(i))^2 + (sqrt(sum(integral_sqerror))/sum(integral_bs))^2)^0.5;
end


CS_fit = [CS(:,1) CS_bs(:,2)]
LW_fit = [LW(:,1) LW_bs(:,2)]
integral_fit = [integral integral_bs_err]
populations_fit = [populations_frac populations_frac_bs_err]


results_fit_total = [CS_fit LW_fit integral_fit./1000 populations_fit]



