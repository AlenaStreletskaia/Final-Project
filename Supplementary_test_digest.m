%% Functions used for the restriction analysis:
%For all bands (that are real to separate on the ~1-2% agarose gel: ~0.4-6kb)
func = @(x) -(20/(0.45*x + 0.75)^8 - 20/(0.8*x+0.5)^2);
fplot(func, [0.4 6], 'LineWidth', 2); 
xlabel('Distance, kb'); ylabel('Separation optimum');
yline(0.94, '-.b'); hold off;
legend('Values of the function', 'Threshold');

%% Maximal value:
x = linspace(0, 5, 100);
func = @(x) -(20/(0.45*x + 0.75)^8 - 20/(0.8*x+0.5)^2);
y = zeros(1,length(x));
for i = 1:length(x)
    y(i) = func(x(i));
end
max = max(y)
%max =
%    7.1899
thres = 0.94/max*100 %threshold in pct from max value
%thres =
%   13.0739