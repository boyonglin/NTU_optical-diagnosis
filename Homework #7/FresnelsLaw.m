function [r, t] = FresnelsLaw(cz, np, nc)
    theta_1 = acos(-cz);                %incident angle
    theta_t = asin(sin(theta_1)*np/nc); %refraction angle
    apha = theta_1 - theta_t;
    beta = theta_1 + theta_t;
    theta_c = asin(nc/np);  %critical angle
    if theta_1 >= theta_c   %total internal reflection
        r = 1; t = 0;
    else
        r = 0.5*(sin(apha)^2/sin(beta)^2 + tan(apha)^2/tan(beta)^2);    %Fresnel’s law
        t = 1-r;    %transmittance
    end
end