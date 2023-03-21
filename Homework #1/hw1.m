% Bo-Yong Lin, 2022/9/7, Random Number Generation

clear;  % remove items from workspace, freeing up system memory


N = 10000;      % number of photons
r = rand(N,1);  % rand generates a uniform random number distributed between (0~1)

t = tiledlayout('flow');            % create a tiled chart layout that can accommodate any number of axes
nexttile                            % create the first axes
h1 = histogram(rand(N,1), 20);      % bar graph simulation result of r.n.
title('Run 1')                      % create a subtitle in each tile

nexttile
h2 = histogram(rand(N,1), 20);
title('Run 2')

nexttile
h3 = histogram(rand(N,1), 20);
title('Run 3')

nexttile
h4 = histogram(rand(N,1), 20);
title('Run 4')

nexttile
h5 = histogram(rand(N,1), 20);
title('Run 5')

% display a shared title and axis labels
title(t, 'Generate random numbers')
xlabel(t, 'Intervals')
ylabel(t, '# of photons')

vect = [h1.Values; h2.Values ; h3.Values ; h4.Values ; h5.Values];
mean = mean(vect);  % calculate mean value of array
std = std(vect);    % calculate standard deviation of array
interval = linspace(0.05, 1, 20);

format bank         % set output display format 2 digits after the decimal point
T = table(interval', mean', std', 'VariableNames', {'Intervals' 'Mean' 'Standard deviation'});
disp(T);            % disply table array with named variables