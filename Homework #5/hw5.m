% Bo-Yong Lin, 2022/10/27, develop a grid structure with infinitely narrow beam

clear; close;

simulationRun = 5;  % number of simulation runs
photonNum = 10000;  % number of photons per run

depth = 0.15;       % tissue depth [cm]
width = 0.3;        % tissue width [cm]
dz = 0.01;          % depth resolution [cm]
dw = 0.01;          % width resolution [cm]

us = 414;           % scattering coefficient [1/cm]
ua = 6;             % absorption coefficient [1/cm]
ut = ua + us;       % extinction coefficient [1/cm]
g = 0.91;           % anisotropy factor
n_air = 1.0;        % refractive index of air
n_tissue = 1.37;    % refractive index of tissue

threshold = 0.1;    % if the weight is significantly small
rouletteIndex = 20; % re-energize photon

Rt = zeros(simulationRun, 0);   % initialize the total reflected array
Tt = zeros(simulationRun, 0);   % initialize the total transmitted array
At = zeros(simulationRun, 0);   % initialize the total absorbed array

for run = 1:simulationRun
    Rsp = 0; Rd = 0; T = 0; A = 0;
    Absorption = zeros(width/dw, depth/dz);  % build Absorption matrix
    for N = 1:photonNum

        % initiallize photon
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;
        terminate = false;
        randNum = rand(1);
        stepSize = -log(randNum)/ut;
        weight = 1.0;

        % first partial reflection
        Rsp = (n_air - n_tissue)^2/(n_air + n_tissue)^2;    % specular reflectance
        weight = (1 - Rsp)*weight;

        % first photon-tissue interaction
        x = x + cx*stepSize;
        y = y + cy*stepSize;
        z = z + cz*stepSize;
        A = A + weight*ua/ut;   % separated absorption due to the first photon-tissue interaction
        weight = weight*us/ut;
        [stepSize, cx, cy, cz] = NewDirection(ut, g, cx, cy, cz);

        while terminate == false

            % update the new position of the photon
            x = x + cx*stepSize;
            y = y + cy*stepSize;
            z = z + cz*stepSize;

            if z > depth
                [r, t] = FresnelsLaw(cz, n_tissue, n_air);

                T = T + t*weight;       % update transmittance
                weight = r*weight;
                z = 2*depth - z;
                cz = -cz;

            elseif z < 0
                [r, t] = FresnelsLaw(-cz, n_tissue, n_air);

                Rd = Rd + t*weight;     % update diffuse reflectance
                weight = r*weight;
                z = -z;
                cz = -cz;
            end

            radius = sqrt(x^2 + y^2);           % current radius
            ir = 1 + round(radius/dw - 0.5);    % grid of width
            iz = 1 + round(z/dz - 0.5);         % grid of depth
            
            % grid bookkeeping
            if (ir <= 30) && (iz <= 15)
                Absorption(ir, iz) = Absorption(ir, iz) + weight*ua/ut;
                A = A + weight*ua/ut;
                weight = weight*us/ut;

                Volume(ir) = (2*ir - 1)*pi*dz*dw^2;                             % the volume of the Ith element
                Source(ir, iz) = Absorption(ir, iz)/(Volume(ir)*photonNum);     % the source term for the current grid element
                Fluence(ir, iz) = Source(ir, iz)/ua;                            % the fluence rate for the grid element
            end

            [weight, terminate, stepSize, cx, cy, cz] = roulette(g, weight, threshold, rouletteIndex, terminate, stepSize, ut, cx, cy, cz);
        end
    end

    Rt(run) = Rsp + Rd/(photonNum - 1);
    Tt(run) = T/photonNum;
    At(run) = A/photonNum;
end

rowNames = ["Run " + num2cell(1:5) "Mean" "Std"];
varNames = {'R' 'T' 'A'};
Rt = [Rt, mean(Rt), std(Rt)];
Tt = [Tt, mean(Tt), std(Tt)];
At = [At, mean(At), std(At)];

format default
T = table(Rt', Tt', At', 'VariableNames', varNames, 'RowNames',rowNames);
disp(T);

tiledlayout('flow')
b1 = nexttile;
bar3(Absorption)
title('Photons Absorbed within the Tissue')
xlabel('Depth [cm]')
ylabel('Width [cm]')
zlabel('# of Photons Absorbed')
xlim([0 15])
ylim([0 30])
xticks([0 15])
yticks([0 30])
xticklabels({'0', '0.15'})
yticklabels({'0', '0.3'})
b1.XDir = 'reverse';

b2 = nexttile;
bar3(Fluence)
title('Fluence Rate of Scarttered Photons')
xlabel('Depth [cm]')
ylabel('Width [cm]')
zlabel('Fluence rate [1/cm^{2}]')
linkaxes([b1 b2], 'xy')
xticks([0 15])
yticks([0 30])
xticklabels({'0', '0.15'})
yticklabels({'0', '0.3'})
b2.XDir = 'reverse';

%% Get new direction for current photon

function [stepSize, cx, cy, cz] = NewDirection(ut, g, cx, cy, cz)
    randNum = rand(1);
    stepSize = -log(randNum)/ut;
    randNum = rand(1);
    phi = 2*pi*randNum;
    randNum = rand(1);
    temp = (1 - g*g)/(1 - g + 2*g*randNum);
    theta = acos(0.5*(1 + g*g - (temp^2))/g);
    
    old_cx = cx; old_cy = cy; old_cz = cz;
    
    % update the directional cosines
    if abs(old_cz) > 0.99999
        cx = sin(theta) * cos(phi);
        cy = sin(theta) * sin(phi);
        cz = old_cz * cos(theta) / abs(old_cz);
    else
        cx = sin(theta)/sqrt(1 - old_cz*old_cz)*(old_cx*old_cz*cos(phi) - old_cy*sin(phi)) + old_cx*cos(theta);
        cy = sin(theta)/sqrt(1 - old_cz*old_cz)*(old_cy*old_cz*cos(phi) + old_cx*sin(phi)) + old_cy*cos(theta);
        cz = -sin(theta)*cos(phi)*sqrt(1 - old_cz*old_cz) + old_cz*cos(theta);
    end
end

%% Calculate partial reflectance (r) and transmittance (t)

function [r, t] = FresnelsLaw(cz, n_in, n_out)
    theta_1 = acos(cz);                         % incident angle
    theta_t = asin(sin(theta_1)*n_in/n_out);    % refraction angle
    theta_c = asin(n_out/n_in);                 % critical angle for totally internally reflected

    % is total internal reflection
    apha = theta_1 - theta_t;
    beta = theta_1 + theta_t;
    
    if theta_1 >= theta_c
        r = 1; t = 0;
    else
        r = 0.5*(sin(apha)^2/sin(beta)^2 + tan(apha)^2/tan(beta)^2);    %Fresnel’s law
        t = 1 - r;
    end
end

%% Let's play roulette!

function [weight, terminate, stepSize, cx, cy, cz] = roulette(g, weight, threshold, rouletteIndex, terminate, stepSize, ut, cx, cy, cz)
    if weight < threshold
        rn = rand(1);
        if rn <= 1/rouletteIndex
            weight = weight*rouletteIndex;
        else
            weight = 0;
            terminate = true;
        end
    else
        [stepSize, cx, cy, cz] = NewDirection(ut, g, cx, cy, cz);
    end
end
