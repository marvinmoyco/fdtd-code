%========== FDTD1D.m, Dec 2007 ============
eps0 = 8.82e-12;
mu0 = pi*4e-7;
eta = sqrt(mu0/eps0);
c = 1./sqrt(eps0*mu0);
Exoldr = 0;
Exoldl = 0;
rbc = 0; %Reflecting boundary at right
lbc = 1; %Reflecting Boundary at left

nt = 500; %input('Number of time steps = ?');
zmax = 1; %input('Maximum distance [m] = ?');
ke = 101; %input('Number of range steps, kc [] = ?');
kc = 30;  %input('Source location, kc [] = ?');

delz = floor(zmax/ke);
dt = delz/c;
tt = ke/2;
t0 = 10;
ce = dt/(1*delz);
ch = dt/(1*delz);

band = c/(10*delz);
alfa = 3.3*band*band;
shift = 4./sqrt(alfa);

Ex = zeros(ke);
Hy = zeros(ke);
z = zeros(ke);


%(1) put zeros (2) Inject a plane wave
for k=1:ke %Initial values along z
    Ex(k) = 0.0;
    Hy(k) = 0.0;
    z(k) = k; %(1)
    %Ex(k) = exp(-(kc-k)^2/tt);
    %Hy(k) = exp(-(kc-k)^2/tt)/eta; %(2)
end

t = 0;
for n=1:nt %FDTD time loop
    t=n*dt;
    for k=2:ke-1 %Ex is calculated along z
        Ex(k) = Ex(k) + ce*(Hy(k-1) - Hy(k));
    end
    
    %Inject Gaussian Pulse
    
    Ex(kc) = Ex(kc) + exp(-alfa*(t - shift)^2);
    
    if lbc == 0
        Ex(1) = Exoldl;
    end
    if rbc==0
        Ex(ke) = Exoldr;
    end
    
    
    for k=1:ke-1 %Hy is calculated along z
        Hy(k) = Hy(k) + ch*(Ex(k) - Ex(k+1));
    end
    Exoldr = Ex(ke-1);
    Exoldl = Ex(2);
    
    plot(z,Ex);
    xlim([1,ke]); 
    ylim([-1.,1.]);
    xlabel('z axis');
    grid;
    pause(0.001);
end %end of loop








