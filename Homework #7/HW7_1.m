clear; close;

runs = 5;   %times of simulated
Ns = 10000; %number of photons per run

d1 = 0.005;     %layer1 depth[cm]
d2 = 0.205;     %layer2 depth[cm]
w  = 0.3;       %tissue width[cm]
dz = 0.0025;    %depth resolution[cm]
dr = 0.0025;    %width resolution[cm]

%refractive index of air/layer1/layer2
n0 = 1.0;
n1 = 1.4;
n2 = 1.4;

%layer1 coefficient[cm^(-1)]
ua1 = 37;
us1 = 480;
ut1 = ua1+us1;
g1 = 0.79;

%layer2 coefficient[cm^(-1)]
ua2 = 2.2;
us2 = 220;
ut2 = ua2+us2;
g = 0.79;

%if the weight is significantly small, re-energize every 1/m photons
threshold = 0.1;
m = 20;

for run = 1:runs
    T = 0; Rsp = 0; Rd = 0; A = 0;
    fc = zeros(1,4); fp = zeros(1,4); %current factor / previous factor [n, ut, ua, us]
    Ab = zeros(w/dr, d2/dz);  %build Absorption matrix
    for N = 1:Ns
        
        %initiallize photon
        x = 0; y = 0; z = 0;
        cx = 0; cy = 0; cz = 1;
        terminate = false;
        weight = 1.0;
        
        %first partial reflection
        Rsp = (n0 - n1)^2/(n0 + n1)^2;
        weight = (1 - Rsp)*weight;
        
        fc = [n1, ut1, ua1, us1];
        
        %first photon-tissue interaction
        rn = rand(1);
        s = -log(rn)/fc(2);
        z = z + cz*s;
        A = A + weight*fc(3)/fc(2);
        weight = weight*fc(4)/fc(2);
        [cx, cy, cz] = NewDirection(g, cx, cy, cz);
        
        while terminate == false
            prez = z;   %previous z position
            rn = rand(1);
            s = -log(rn)/fc(2);
            %update the new position of the photon
            x = x + cx*s;
            y = y + cy*s;
            z = z + cz*s;
            
            n = 0; %n means no
            while n < 2
                if d2 < z
                    [r, t] = FresnelsLaw(-cz, n2, n0);
                    T = T + t*weight;   %transmittance
                    weight = r*weight;
                    z = 2*d2-z;
                    cz = -cz;

                elseif z < 0
                    [r, t] = FresnelsLaw(cz, n1, n0);
                    Rd = Rd + t*weight; %diffuse reflectance
                    weight = r*weight;
                    z = -z;
                    cz = -cz;
                else
                    n = 1;
                end

                if (0 <= prez) && (prez < d1) && (d1 <= z) && (z < d2) %into another layer? yes
                    prez = z;
                    x = x - cx*s;
                    y = y - cy*s;
                    z = z - cz*s;

                    s1 = (d1-z)/cz;
                    s2 = ut1/ut2*(s-s1);
                    x = x + cx*(s1+s2);
                    y = y + cy*(s1+s2);
                    z = z + cz*(s1+s2);
                elseif (d1 <= prez) && (prez < d2) && (0 <= z) && (z < d1)
                    prez = z;
                    x = x - cx*s;
                    y = y - cy*s;
                    z = z - cz*s;

                    s2 = (d1-z)/cz;
                    s1 = ut2/ut1*(s-s2);
                    x = x + cx*(s1+s2);
                    y = y + cy*(s1+s2);
                    z = z + cz*(s1+s2);
                else
                    n=2;
                end
            end
            
            if z < 0
                fc = [n0, ut1, ua1, us1];
            elseif (0 <= z) && (z < d1)
                fc = [n1, ut1, ua1, us1];
            elseif (d1 <= z) && (z < d2)
                fc = [n2, ut2, ua2, us2];
            elseif d2 <= z
                fc = [n0, ut2, ua2, us2];
            end
            
            
            ra = sqrt(x^2 + y^2);        %current radius
            ir = 1+round(ra/dr - 0.5);   %grid of width
            iz = 1+round(z/dz - 0.5);    %grid of depth
            
            if  (ir<=w/dr) && (iz<=d2/dz)
                Ab(ir,iz) = Ab(ir,iz) + weight*fc(3)/fc(2);
                A = A + weight*fc(3)/fc(2);
                weight = weight*fc(4)/fc(2);

                V(ir) = (2*ir-1)*pi*dz*dr^2;        %the volume of the Ith element 
                S(ir,iz) = Ab(ir,iz)/(V(ir)*Ns);    %the source term for the current grid element 
                F(ir,iz) = S(ir,iz)/fc(3);          %the fluence rate for the grid element
            end
            
            if weight < threshold   %roulette
                rn = rand(1);
                if rn <= 1/m
                    weight = weight*m;
                    [cx, cy, cz] = NewDirection(g, cx, cy, cz);
                else
                    weight = 0;
                    terminate = true;
                end
            else
                [cx, cy, cz] = NewDirection(g, cx, cy, cz);
            end
        end
    end

    Rt = Rsp + Rd/Ns;
    Tt = T/Ns;
    At = A/Ns;
    Ra(run) = Rt; Ta(run) = Tt;  %array used to accumulate
    fprintf('\n');
    fprintf('Run %d　　R = %.4f　　T = %.4f　　A = %.4f', run, Rt, Tt, At);
end

subplot(1,2,1);
bar3(Ab);
xlabel('Depth')
ylabel('Width')
zlabel('# of Photons Absorbed')
title('Photons Absorbed within the Tissue: Run(5/5)')
xlim([0 82])
ylim([0 40])
xticks(0:10:82)
yticks(0:10:40)
ax = gca;
ax.XDir = 'reverse';


subplot(1,2,2);
bar3(F);
xlabel('Depth')
ylabel('Width')
zlabel('Fluence rate [1/cm^{2}]')
title('Fluence Rate of Scarttered Photons: Run(5/5)')
xlim([0 82])
ylim([0 40])
xticks(0:10:82)
yticks(0:10:40)
ax = gca;
ax.XDir = 'reverse';

fprintf('\n\n');
RM = mean(Ra);  TM = mean(Ta);  %means of R, T
RS = std(Ra);   TS = std(Ta);   %standard deviations of R, T
fprintf('=> R = %.4f ± %.4f\n', RM, RS);
fprintf('   T = %.4f ± %.4f\n', TM, TS);
fprintf('\n');