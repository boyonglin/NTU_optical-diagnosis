% Bo-Yong Lin, 2022/10/14, compare the effect of index-matching at difference surfaces

clear;

n_tissue = 1.4;     % refractive index of tissue
n_air = 1.0;        % refractive index of air
n_water = 1.33;     % refractive index of water

medium = {'air' 'water'};
a = ['air' 'water'];
Rd = [];

[first_term, second_term] = rd(n_tissue, n_air);
Rd(1) = first_term + second_term;

[first_term, second_term] = rd(n_tissue, n_water);
Rd(2) = first_term + second_term;

T = table(string(medium)', round(Rd, 4, "significant")', 'VariableNames', {'Medium' 'Rd'});
disp(T);

%% reflection coefficient of diffuse light

function [first_term, second_term] = rd(n_in, n_out)
    theta_c = asin(n_out/n_in);
    f = @(theta) 0.5*(sin(theta - asin(sin(theta)*n_in/n_out)).^2./sin(theta + asin(sin(theta)*n_in/n_out)).^2 + tan(theta - asin(sin(theta)*n_in/n_out)).^2./tan(theta + asin(sin(theta)*n_in/n_out)).^2).*sin(2*theta);
    first_term = integral(f, 0, theta_c);
    f = @(theta) sin(2.*theta);
    second_term = integral(f, theta_c, pi/2);
end