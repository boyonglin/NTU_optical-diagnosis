function [cx, cy, cz] = NewDirection(g, cx, cy, cz)
    rn = rand(1);
    phi = 2*pi*rn;
    rn = rand(1);
    temp = (1-g*g)/(1-g+2*g*rn);
    theta = acos(0.5*(1+g*g-(temp^2))/g);

    old_cx = cx; old_cy = cy; old_cz = cz;

    %update the directional cosines
    if abs(old_cz) > 0.99999
        cx = sin(theta)*cos(phi);
        cy = sin(theta)*sin(phi);
        cz = old_cz*cos(theta)/abs(old_cz);
    else
        cx = sin(theta)/sqrt(1-old_cz*old_cz)*(old_cx*old_cz*cos(phi)-old_cy*sin(phi)) + old_cx*cos(theta);
        cy = sin(theta)/sqrt(1-old_cz*old_cz)*(old_cy*old_cz*cos(phi)+old_cx*sin(phi)) + old_cy*cos(theta);
        cz = -sin(theta)*cos(phi)*sqrt(1-old_cz*old_cz) + old_cz*cos(theta);
    end
end