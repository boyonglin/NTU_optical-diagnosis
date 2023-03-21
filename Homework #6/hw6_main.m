% Bo-Yong Lin, 2022/11/02, plot the fluence rate distribution with finite beam

clear; close;

[fluence_flatTop] = hw6_functions('flatTop');
[fluence_gaussian] = hw6_functions('gaussian');

tiledlayout('flow')
b1 = nexttile;
bar3(fluence_flatTop)
title('Fluence Rate of Scarttered Photons with Uniform flat-top beam')
xlabel('Depth [cm]')
ylabel('Width [cm]')
zlabel('Fluence rate [1/cm^{2}]')
xlim([0 15])
ylim([0 30])
xticks([0 15])
yticks([0 30])
xticklabels({'0', '0.15'})
yticklabels({'0', '0.3'})
b1.XDir = 'reverse';

b2 = nexttile;
bar3(fluence_gaussian)
title('Fluence Rate of Scarttered Photons with Gaussian beam')
xlabel('Depth [cm]')
ylabel('Width [cm]')
zlabel('Fluence rate [1/cm^{2}]')
linkaxes([b1 b2], 'xy')
xticks([0 15])
yticks([0 30])
xticklabels({'0', '0.15'})
yticklabels({'0', '0.3'})
b2.XDir = 'reverse';