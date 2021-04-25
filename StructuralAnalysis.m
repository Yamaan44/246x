
%% Inputs
clear; close

%World constants
g = 9.81; %m/s^2, gravitational constant
rho = 1.225; %kg/m^3, air density

%Wing Geometry
AR = 5; %aspect ratio
S = 0.96; %m^2, wing surface area
lambda = 0.45; %taper ratio

rho_f = 25.23; %oz/cu-ft, Foam density, potentially 40 
rho_f = rho_f* 1.00115; %oz/cu-ft to kg/m^3

%Horizontal Tail Geometry
AR_t = 3; %aspect ratio
S_t = 0.3; %m^2, horizontal stab surface area
lambda_t = 0.45; %taper ratio

rho_f_t = 25.23; %oz/cu-ft, Foam density, potentially 40 
rho_f_t = rho_f_t* 1.00115; %oz/cu-ft to kg/m^3

%Vertical Tail Geometry
AR_tv = 2; %aspect ratio
S_tv = 0.15; %m^2, vertical stab surface area
lambda_tv = 0.45; %taper ratio

rho_f_tv = 25.23; %oz/cu-ft, Foam density, potentially 40 
rho_f_tv = rho_f_tv* 1.00115; %oz/cu-ft to kg/m^3

%Spar properties
E = 7200; %MPa 
I = 169412 * 0.001^4; %m^4, derived from spar geometry
J = 169412 * 0.001^4; %m^4, derived from spar geometry
ymax = 0.02; %m, distance from neutral axis
r = 0.02; %m, torsional "radius" from centroid
d = 0.01; %m center of pressure to spar distance

%Boom/Fuselage properties
E_b = 7200; %MPa
I_b = 169412 * 0.001^4; %m^4, derived from spar geometry
J_b = 169412 * 0.001^4; %m^4, derived from spar geometry
ymax_b = 0.02; %m, distance from neutral axis
r_b = 0.02; %m, torsional "radius" from centroid
rho_boom = 0.046; %kg/m, 10cm x 10cm square

%Aircraft Geometry
s = 0.3; %m, empennage cp to wing cp
s_cg = 0.05; %m, cg to wing cp (cg behind cp)

%Forces
L_S = 25; %N/m^2, wing loading
L_min = 2; %N, turning
L_max = 10;%N, turning
L_e = 1; %N, empennage lift

%Limits
d_max = 0.02; %m, max tip deflection

%masses from spreadsheet
m_payload = 0.7513; %kg


%% Derived Inputs
b = sqrt(S*AR); %m, span
m_boom = rho_boom * b; %kg, constant
cbar = sqrt(S/AR); %m, mean chord

cr = 2*S/(b*(1+lambda));
ct = cr * lambda;

bt = sqrt(S_t*AR_t);
crt = 2*S_t/(bt*(1+lambda_t));
ctt = crt*lambda_t;

btv = sqrt(S_tv*AR_tv);
crtv = 2*S_tv/(btv*(1+lambda_tv));
cttv = crtv*lambda_tv;

xvec = linspace(0,b/2,10000); %m, x positions along wing starting from root to b/2
cvec = ((b/2 - xvec)*cr + xvec*ct)*2/b; %m, chord at each x position along wing starting from root to b/2

L = L_S*S; %N, total lift

%%% USE AIRFOIL READER TO CALCULATE WING VOLUME, THEN MASS
[~,xin,yin]=openfile('naca4412.dat');
plot(xin,yin);
a = polyarea(xin,yin); %m^2/m (chord)
vol = 2*trapz(xvec,a*cvec); %m^3, total wing volume
m_wing = vol * rho_f + rho_boom * b; %kg, wing mass (assumes solid x section)



%TAIL WEIGHT, HORIZONTAL
xtvec = linspace(0,bt/2,10000); %m, x positions along tail starting from root to bt/2
ctvec = ((bt/2 - xtvec)*crt + xtvec*ctt)*2/bt; %m, chord at each x position along tail starting from root to b/2

[~,xtin,ytin]=openfile('naca4412.dat');
plot(xtin,ytin);
at = polyarea(xtin,ytin); %m^2/m (chord)
volt = 2*trapz(xtvec,at*ctvec); %m^3, total tail volume

%TAIL WEIGHT, VERTICAL
xtvvec = linspace(0,btv,10000); %m, x positions along tail starting from root to btv
ctvvec = ((btv - xtvvec)*crtv + xtvvec*cttv)/btv; %m, chord at each x position along tail starting from root to b/2

[~,xtin,ytin]=openfile('naca4412.dat');
plot(xtin,ytin);
atv = polyarea(xtin,ytin); %m^2/m (chord)
voltv = trapz(xtvvec,atv*ctvvec); %m^3, total tail volume

m_tail = volt * rho_f_t + voltv * rho_f_tv + btv* rho_boom + bt * rho_boom; %kg, tail horizontal mass (assumes solid x section)

%% Mass

m_total = m_wing + m_tail + m_boom + m_payload;
%% Tail boom forces
W = m_total*g;
x_bvec = linspace(0,s,10000); %m, boom coord sys, origin at cp of wing
V_boom = zeros(10000,1);
idx = find(x_bvec>s_cg,1);
V_boom(1:idx) = L;
V_boom(idx+1:end) = L-W;
M_boom = cumtrapz(x_bvec,V_boom);


u_primevec_b = cumtrapz(x_bvec,M_boom);
u_bvec = E*I*cumtrapz(x_bvec,u_primevec_b); %m
u_bmax = max(u_bvec);

figure(4)
subplot(3,1,1)
plot(x_bvec,V_boom);
xlabel('x (m)')
ylabel('V(x) (N)')
subplot(3,1,2)
plot(x_bvec,M_boom);
xlabel('x (m)')
ylabel('M(x) (Nm)')
subplot(3,1,3)
plot(x_bvec,u_bvec);
xlabel('x (m)')
ylabel('u(x) (m)')

%% Deflection Calculation

%chord as a fcn of span
figure(1)
fplot(@(x) ((b/2 - x)*cr + x*ct)*2/b, [0 b/2])

xlim([0 b/2]);
ylim([0 cr]);

%wvec = cvec * L_S; %N/m, rectangular lift distribution, do not use
w0 = 4*L/(pi*b);
wvec = w0 * sqrt(1 - (xvec/(0.5*b)).^2);%N/m, elliptical lift distribution


Vvec = cumtrapz(xvec,wvec); %N
Vvec = flip(Vvec);

Mvec = cumtrapz(xvec,Vvec); %Nm
Mvec = flip(Mvec);

u_primevec = cumtrapz(xvec,Mvec);
uvec = E*I*cumtrapz(xvec,u_primevec); %m
umax = max(uvec);


figure(2)
fplot(@(x) ((b/2 - x)*cr + x*ct)*2/b * L_S, [0 b/2])

figure(3)
subplot(3,1,1)
plot(xvec,Vvec);
xlabel('x (m)')
ylabel('V(x) (N)')
subplot(3,1,2)
plot(xvec,Mvec);
xlabel('x (m)')
ylabel('M(x) (Nm)')
subplot(3,1,3)
plot(xvec,uvec);
xlabel('x (m)')
ylabel('u(x) (m)')

%% Wing Stresses

V_max = max(Vvec);
M_max = max(Mvec);
T_spar = L * d;

sigma_max = M_max * ymax/I * 10^-6; %MPa, max compressive/tensile at root
tau_max =T_spar * r/J * 10^-6; %MPa, max shear due to torsion

%% Boom stresses

T_boom = (L_max - L_min)*0.2; %Nm, torque due to turn %more fidelity into center of lift needed

sigma_max_boom = max(M_boom) * ymax_b/I_b * 10^-6; %MPa, max compressive/tensile at root

tau_max_boom =T_boom * r_b/J_b * 10^-6; %MPa, max shear due to torsion

%% Outputs

fprintf('Wing Tip deflection is %f m \n', umax)
fprintf('Max Wing tensile/compressive stress is %f MPa\n', sigma_max)
fprintf('Max Boom shear stress is %f MPa\n', tau_max)
fprintf('Max Boom deflection is %f m \n', u_bmax)
fprintf('Max Boom tensile/compressive stress is %f MPa\n', sigma_max_boom)
fprintf('Max Boom shear stress is %f MPa\n', tau_max_boom)
fprintf('Total mass is %f kg\n', m_total)
%fprintf('Max shear stress is %f', sigma_max)




function [n,xin,yin] = openfile(file)
    % Read camber line x and y coordinate data from file ?camber.dat?
    % and place in vectors xin and yin
    fid = fopen(file);
    n = fscanf(fid,'%10f',1); %length of airfoil vectors
    data = fscanf(fid,'%10f %10f \n', [2,n]);
    fclose(fid);
    xin = data(1,:);
    big=max(xin);
    xin=xin./big;
    yin = data(2,:);
    yin=yin./big;
end