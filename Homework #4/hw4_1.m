% Bo-Yong Lin, 2022/10/5, Index Mismatch of Boundary Layers (variable-weight)

clear;

d = 0.2;            % tissue depths [cm]
us = 90;            % scattering coefficient [cm^-1]
ua = 10;            % absorption coefficient [cm^-1]
ut = ua + us;       % extinction coefficient
n0 = 1.0;           % refractive index of air
n1 = 1.5;           % refractive index of tissue

threshold = 0.1;    % if the weight is significantly small
m = 20;             % re-energize every 1/m photons

run = 5;            % number of simulations
Rt = zeros(run,0);  % initialize the Rt array

for runs = 1:run
    T = 0; Rsp = 0; Rd = 0; A = 0;  % Transmittance; SPecular Reflection; Diffuse Reflectance; Absorption
    for N = 1:10000

        % initiallize photon
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;
        terminate = false;
        rn = rand(1);
        s = -log(rn)/ut;
        weight = 1.0;

        % first partial reflection
        Rsp = (n0 - n1)^2/(n0 + n1)^2;    % specular reflectance, when the beam is normal to the surface
        weight = (1 - Rsp)*weight;        % update photon weight

        while terminate == false

            % update the new position of the photon
            x = x + cx*s;
            y = y + cy*s;
            z = z + cz*s;

            if z > d
                theta_i = acos(cz);                 % incident angle
                theta_t = asin(sin(theta_i)*n1/n0); % refraction angle
                theta_c = asin(n0/n1);              % critical angle for totally internally reflected

                [r, t] = TIR(theta_i, theta_t, theta_c);

                T = T + t*weight;   % update transmittance
                weight = r*weight;  % update weight
                z = 2*d - z;        % update z coordinate
                cz = -cz;           % update cz

                % absorption of photons reflected from the surface back to the tissue
                A = A + weight*ua/ut;
                weight = weight*us/ut;


            elseif z < 0
                theta_i = acos(-cz);
                theta_t = asin(sin(theta_i)*n1/n0);
                theta_c = asin(n0/n1);

                [r, t] = TIR(theta_i, theta_t, theta_c);

                Rd = Rd + t*weight; % update diffuse reflectance
                weight = r*weight;  % update weight
                z = -z;             % update z coordinate
                cz = -cz;           % update cz

                %absorption of photons reflected from the surface back to the tissue
                A = A + weight*ua/ut;
                weight = weight*us/ut;

            else
                A = A + weight*ua/ut;
                weight = weight*us/ut;
            end

            [weight, terminate, s, cx, cy, cz] = roulette(weight, threshold, m, terminate, s, ut, cx, cy, cz);

        end
    end

    Rt(runs) = Rsp + Rd/9999;   % total reflectance
end

runs = 1:run;
format default
T = table(runs', Rt', 'VariableNames', {'Run' 'Rt'});
disp(T);

%% is total internal reflection

function [r, t] = TIR(theta_1, theta_t, theta_c)
    apha = theta_1 - theta_t;
    beta = theta_1 + theta_t;

    if theta_1 >= theta_c
        r = 1; t = 0;
    else
        r = 0.5*(sin(apha)^2/sin(beta)^2 + tan(apha)^2/tan(beta)^2);    % Fresnel's law
        t = 1 - r;      % partial transmittance
    end
end


%% let's play roulette!

function [weight, terminate, s, cx, cy, cz] = roulette(weight, threshold, m, terminate, s, ut, cx, cy, cz)
    if weight < threshold
        rn = rand(1);
        if rn <= 1/m
            weight = weight*m;
        else
            weight = 0;
            terminate = true;
        end
    else
        % calculate new direction
        rn = rand(1);
        s = -log(rn)/ut;
        rn = rand(1);
        phi = 2*pi*rn;
        rn = rand(1);
        theta = acos(2*rn - 1);
    
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
end