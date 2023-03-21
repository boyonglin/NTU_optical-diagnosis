% Bo-Yong Lin, 2022/9/14, Attenuation of a Collimated Beam

clear;  % remove items from workspace, freeing up system memory

N = 10000;      % number of photons
ua = 10;        % absorption coefficient [cm^-1]
dz = 0.025;     % depth interval [cm]

depthInterval = 0:dz:1;                 % depth interval as x-axis
tiledChart = tiledlayout('flow');       % create a tiled chart layout that can accommodate any number of axes

for times = 1:5
    nexttile                                    % create the first axes
    stepsize=(-log(rand(N, 1)))/ua;             % variable step method
    histogram(stepsize, depthInterval);         % bar graph simulation result of r.n.
    xlim([0, 1]);                               % set x-axis lmits
    xticks(0:0.1:1);                            % set x-axis tick values
    hold on
    f = 10000*dz*ua*exp(-ua*depthInterval);     % the solid line is that from Beer-Lambert law
    plot(depthInterval, f)
    title1 = "Run ";
    title2 = num2str(times);
    title(strcat(title1, title2))       % concatenate strings horizontally
end

% display a shared title and axis labels
title(tiledChart, 'Attenuation of a Collimated Beam')
xlabel(tiledChart, 'Intervals')
ylabel(tiledChart, '# of photons absorbed')