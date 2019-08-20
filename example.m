clear
data = randn(10000, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set optional parameters.  
% Default values listed in brackets []
% Allowed values listed in braces {}
% To accept all defaults, call EstimatePDF without parameters argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters.SURDtarget = 70;                 % [40]          {number > 0 and <= 100}
parameters.SURDmin = 50;                    % [5]           {number > 0}
parameters.SURDmax = 85;                    % [100]         {number > SURDmin and <= 100}
parameters.LagrangeMin = 3;                 % [1]           {number >= 1}
parameters.LagrangeMax = 10;                % [200]         {number > LagrangeMin}
parameters.lowBound = -5;                   % [-infinite]   {any number}
parameters.highBound = 5;                   % [infinite]    {number > lowBound}
parameters.integrationPoints = 200;         % [calculated]
parameters.partition = 0;                   % [1024]        {any number, zero for no partitioning}
parameters.fuzz = false;                    % [false]       {true, false} *untested functionality*
parameters.debug = true;                    % [false]       {true, false}
parameters.scoreType = 'LL';                % ['QZ']        {'QZ', 'LL', 'AD'}
parameters.minVariance = true;              % [true]        {true, false}
parameters.outlierCutoff = 0;               % [7]           {number >= 0, zero to keep all outliers}
parameters.adaptiveDx = true;               % [false]       {true, false}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Call PDF Estimator, catch and display errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    [failed, x, pdf, cdf, sqr, lagrange, score, confidence, SURD] = EstimatePDF(data, parameters);
catch ME
    failed = true;
end

if failed
    disp(ME.message);
    return;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Total number of Lagrange multipliers: ', num2str(length(lagrange))]);
disp(['Threshold score: ', num2str(score)]);
disp(['Confidence level: ', num2str(confidence), '%']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set properties for high quality figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
width = 6;                                                                 % Width in inches
height = 4;                                                                % Height in inches
alw = 1.5;                                                                 % AxesLineWidth 
fsz = 14;                                                                  % Fontsize 
lw = 1.5;                                                                  % LineWidth 
msz = 8;                                                                   % MarkerSize 
set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);   
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw);
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
set(0,'defaultFigurePosition', [400, 50, width*100, height*110]); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Probability Distribution Function (PDF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
plot(x, pdf, 'color', 'black');
xlabel('Sample Data');
ylabel('PDF');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Cumulative Distribution Function (CDF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
plot(x, cdf, 'color', 'black');
xlabel('Sample Data');
ylabel('CDF');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Scaled Quantile Residual (SQR) with lemon drop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smallN = 256;
smallN2 = 258;
rangeSQR = 0:1/(smallN+1):1;       
mu = rangeSQR*(smallN + 1) / (smallN + 1);
lemonDrop = sqrt(mu.*(1-mu)) * 3.4;
sampleCount2 = (smallN + 2):-1:1;
colorRange = (255-204)*sampleCount2/(smallN + 2);
base = repmat(204, smallN + 2, 1);
col = (base + colorRange') / 255;
col255 = repmat(255, smallN + 2, 1) / 255;
rgb = [col255 col255 col];

figure;
hold on;
count2 = 1;
for ii = ceil(smallN2/2):smallN2-1
    ix = [ii ii+1 smallN2-ii smallN2-ii+1];
    h1 = fill(rangeSQR(ix), lemonDrop(ix), rgb(count2, :),'edgecolor', 'none');
    h2 = fill(rangeSQR(ix), -lemonDrop(ix), rgb(count2, :),'edgecolor', 'none');
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    count2 = count2 + 2;
end
xPosition = 0:1/(length(sqr) - 1):1;
plot(xPosition, sqr, 'color', 'black');
xlabel('Position');
ylabel('SQR');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot Sampled Uniform Random Data (SURD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
plot(xPosition, SURD, 'color', 'black');
xlabel('Position');
ylabel('SURD');


