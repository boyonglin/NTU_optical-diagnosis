% Bo-Yong Lin, 2022/10/13, find critical angle from air (n0) into tissue (n1)

clear; close;           % close all open figures

n0 = 1.0;               % refractive index of air
n1 = 1.4;               % refractive index of tissue
row1 = linspace(0, pi/2, 1000);
row2 = zeros(0);
r = [row1; row2];       % initialize the array of specular reflectance
colume = 1;             % counter for r

for theta_i = linspace(0, pi/2, 1000)   % incident angles (0 ~ 90)
    theta_t = asin(sin(theta_i)*n1/n0);
    apha = theta_i - theta_t;
    beta = theta_i + theta_t;
    theta_c = asin(n0/n1);
    if theta_i >= theta_c               % is total internal reflection
        r(2, colume) = 1;
    else                                % be calculated by Fresnelbs law
        r(2, colume) = 0.5*(sin(apha)^2/sin(beta)^2 + tan(apha)^2/tan(beta)^2);
    end
    colume = colume + 1;
end

theta_c_degree = rad2deg(theta_c);      % convert angle from radians to degrees

% mark theta_c
plot(row1, r(2,:), theta_c, 1, 'k*')
text(theta_c, 1, ['\uparrow \theta_c = ' num2str(theta_c_degree) char(176)], 'VerticalAlignment', 'top')

xticks([0 pi/8 pi/4 pi/2])
xticklabels({'0','\pi/8','\pi/4','\pi/2'})
xlabel('incident angle [degrees]')
ylabel('specular reflectance')
title('Specular Reflectance Coefficient')