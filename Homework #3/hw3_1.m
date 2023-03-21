% Bo-Yong Lin, 2022/9/21, Scattering of Photons in Tissue (fixed-weight)

clear;

d = 0.02;       % tissue depths [cm]
us = 90;        % scattering coefficient [cm^-1]
ua = 10;        % absorption coefficient [cm^-1]
ut = ua + us;   % extinction coefficient

for times = 1:5
    T = 0; R = 0; A = 0;    % Transmittance; Reflectance; Absorption
    for N = 1:10000             % number of photons
        x = 0; y = 0; z = 0;    % initial position of photons
        cx = 0; cy = 0; cz = 1; % initial values for the trajectory
        terminate = false;      % determine whether the photon is terminated
        rn = rand(1);
        s = -log(rn)/ut;        % variable step

        while terminate == false
            % update the new position of the photon
            x = x + cx*s;
            y = y + cy*s;
            z = z + cz*s;

            if z > d
                T = T + 1;
                terminate = true;
            elseif z < 0
                R = R + 1;
                terminate = true;
            else
                rn = rand(1);

                if rn <= ua/ut
                    A = A + 1;
                    terminate = true;
                else
                    rn = rand(1);
                    s = -log(rn)/ut;
                    rn = rand(1);
                    phi = 2*pi*rn;          % deflection angle
                    rn = rand(1);
                    theta = acos(2*rn-1);   % azimuthal angle

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
        end
    end
    nexttile
    pie([T R A], '%.2f%%');
    legend('Transmitted', 'Reflected', 'Absorped', 'Location', 'bestoutside')
    titleStr=sprintf('Run %d', times);
    title(titleStr);
end
