function [fitresult, gof] = firndensity_curvefit(density_depths,density_levels,median_firnair)
%CREATEFIT(DENSITY_DEPTHS,DENSITY_LEVELS)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : density_depths
%      Y Output: density_levels
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Apr-2022 10:55:02


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( density_depths, density_levels );

% Set up fittype and options.
ft = fittype( '917-(917-a)*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [350 median_firnair];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'density_levels vs. density_depths', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'density_depths', 'Interpreter', 'none' );
% ylabel( 'density_levels', 'Interpreter', 'none' );
% grid on


