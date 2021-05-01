clc;
close all;
clear all;
% ALL UNITS IN SI %
% Input Parameters to optimize:
% weight - aicraft weight (kg)
% b - wing span (m)
% V_cruise - Cruise velocity (m/s)
global cst % structure with all shared constants
cst = struct();
cst.rho = 1.225; %kg/m^3
cst.g = 9.807; %m/s^2
cst.CL_max = 1.2;
cst.V_stall = 7; % m/s
cst.W_L = 1/2 * cst.rho * cst.V_stall^2 * cst.CL_max / cst.g; % wing loading is sized by stall speed
% Lower and upper bounds
Weight_lo = 0; %Kg
Span_lo = .1; %m
V_cruise_lo = 0; %m/s
Weight_up = 15; %Kg
Span_up = 3; %m
V_cruise_up = 50; %m/s
lb = [Weight_lo,Span_lo,V_cruise_lo];
ub = [Weight_up,Span_up,V_cruise_up];
% Initial conditions
x_0 = [4,1,10];
% Run optimizer
opts = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluation', 10000);
x = fmincon( @(x) optimize( x(1),x(2),x(3)), x_0,[],[],[],[],lb,ub, @(x) nonl_const( x(1),x(2),x(3) ), opts);
% Show optimized aircraft 
weight = x(1)
b = x(2)
v_air = x(3)
AR = b^2/(weight/cst.W_L)
S_ref = weight / cst.W_L;
mass_empty = empty_weight(S_ref, b);
[CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_air);
[mass_bat] = battery_weight(weight, L_D);
[mass_motor,T_W,Power_max] = motor_weight(weight, Drag, v_air);
mass_total = mass_bat + mass_motor + mass_empty;
battery_frac = mass_bat/mass_total
motor_frac = mass_motor/mass_total
structural_frac = mass_empty/mass_total
plot_aircraft(S_ref, b)
%% plot aircraft 
function [] = plot_aircraft(S_ref, b)
  boom_dia = 0.0154; %m - visual purpose only
  [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b);
  %estimate fuselage area
  l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
  r_fus = l_fus./8./2; %assume fuselage fineness ratio of 8
  % we model the fuselage as cone nose and cylinder fuselage
  % assume nose is 1/5 length of entire fuselage
  l_nose = 1/5.*l_fus;
  tail_end = l_t+.75*c_htail;
  pgon = polyshape([tail_end tail_end tail_end-l_fus tail_end-l_fus],[-r_fus r_fus r_fus -r_fus]);
  figure;
  hold on
  rectangle('Position',[(c/4 + l_t-.25*c_htail) (0-(s_htail/c_htail)/2) c_htail s_htail/c_htail],'EdgeColor','r','LineWidth',3); %horizontal tail
  rectangle('Position',[ (c/4 + l_t-.25*c_vtail) 0 c_vtail s_vtail/c_vtail],'EdgeColor','b','LineWidth',3); %vertical tail
  rectangle('Position',[0 0-b/2 c b],'EdgeColor','k','LineWidth',3); %wing
  plot(pgon)
  axis equal
end
%% optimizing functions
function score = optimize(weight, b, v_air)   
  score =-v_air; %maximize velocity given constraints
end
% nonlinear constraints

function [c, ceq] = nonl_const(weight, b, v_air)
  global cst
  S_ref = weight / cst.W_L;
  mass_empty = empty_weight(S_ref, b);
  [CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_air);
  [mass_bat] = battery_weight(weight, L_D);
  [mass_motor,T_W,Power_max] = motor_weight(weight, Drag, v_air);
  mass_total = mass_bat + mass_motor + mass_empty;
  %constrain aspect ratio
  AR = b.^2./S_ref;
  AR_lower = 3.5;
  AR_upper = 20;
%   c = [- AR + AR_lower, bend_stress - max_stress];
  c = [- AR + AR_lower, AR - AR_upper, mass_total - weight];
%   ceq = abs(weight-mass_total);
  ceq = [];
end
%% aircraft sizing functions
function [c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b)
%gives conventional tail dimensions based on wing geometry
%guess some typical values for tail parameters
static_margin = 0.05;
ht_vol_cf = 0.40;
vt_vol_cf = 0.03;
AR_ht = 4;
AR_vt = 1.5; 
c = S_ref./b; %mean chord estimate
l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
l_t = 0.75.*l_fus; %ballpark for length of 1/4 chord to 1/4 tail chord. Should be optimized for weight!
s_htail = ht_vol_cf.*S_ref.*c./l_t;
c_htail = sqrt(s_htail./AR_ht);
s_vtail = vt_vol_cf.*S_ref.*b./l_t;
c_vtail = sqrt(s_vtail./AR_vt);
end 

function [CL, Cd, Cdi , Cd0, L_D, v_ideal, Drag] = calc_aero(weight, b, v_air)
global cst
S_ref = weight / cst.W_L;
[c, s_htail, c_htail, s_vtail, c_vtail, l_t] = size_plane(S_ref, b);
AR = b.^2./S_ref;
CL = 9.81.*cst.W_L./(0.5.*cst.rho.*v_air.^2);
e = .9; %this obviously varies with speed
K = 1./(pi.*AR.*e);
Cdi = CL.^2./(pi.*AR.*e); %lift induced drag
%estimate fuselage area
l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
r_fus = l_fus./8./2; %assume fuselage fineness ratio of 8
%we model the fuselage as cone nose and cylinder fuselage
%assume nose is 1/5 length of entire fuselage
l_nose = 1/5.*l_fus;
Cdf_nose = 2./sqrt(3).*calc_Cf(l_nose,v_air).*(pi.*r_fus.*sqrt(r_fus.^2+l_nose.^2))./S_ref;
Cdf_fuse = calc_Cf(l_fus,v_air).*(2.*pi.*r_fus.*l_fus)./S_ref;
Cdf_wing = calc_Cf(c,v_air).*(2.*S_ref)./S_ref;
Cdf_tail = calc_Cf(c_htail,v_air).*(2*s_htail+2*s_vtail)./S_ref; 
Cd0 = 1.25.*(Cdf_nose + Cdf_fuse + Cdf_wing + Cdf_tail);
Cd = Cd0 + Cdi;
L_D = CL/Cd;
%Find the optimum cruise speed, v_ideal, that maximizes L/D for this
%configuration
Cl_opt = sqrt(Cd0/K);
v_ideal = sqrt( (2*9.81*weight)/(cst.rho*Cl_opt*S_ref) );
Drag = 1/2 * cst.rho * v_air^2 * S_ref * Cd;
  function CF = calc_Cf(l,v_air)
    %flat plate assumption for Cf
    nu = 15.52e-6; %m^2/s
    Re = v_air.*l./nu;
    if Re > 5e5
      CF = 0.455./(log10(Re).^2.58);
    else
      CF = 1.328./sqrt(Re);
    end
  end
end

function [mass_bat] = battery_weight(weight, L_D)
  global cst
  eta_sys = .374; %battery + propellor + motor effiency estimate
  range = 42195/2; %meters
  H_batt = 155*3600; %specific energy density (Joules/kg)
  mass_bat = range.*cst.g.*weight./(L_D.*eta_sys.*H_batt);
end

function [mass_motor,T_W,Thrust,Power_max] = motor_weight(weight, Thrust, v_air)
  Power_cruise = Thrust * v_air; %watts
  Power_max = Power_cruise./0.40; %assume motor runs most efficient at 40% max power
  motor_p_density = 5.101/0.001; %watts/kg of a motor 
  T_W = Thrust*2/(weight*9.81); %static thrust to weight. assume static is twice cruise thrust 
  mass_motor = 1./motor_p_density.*Power_max * 2; % times 2 for prop efficiency
end

function [mass_empty] = empty_weight(S_ref, b)
  [s_htail, c_htail, s_vtail, c_vtail] = size_plane(S_ref, b);
  %estimate fuselage area
  l_fus = 0.75.*b; %ballpark for fuselage length. RC airplanes typically 75% of span
  %weight estimations based off walrus
  m_wing = S_ref.*1.31346914; %kg/m^2 x s_wing
  m_htail = s_htail.*1.695180676; %kg/m^2 x s_htail
  m_vtail = s_vtail.*0.792369773; %kg/m^2 x s_vtail 
  m_fuse = l_fus * 0.3466666667; %kg/m x l_fuselage 
  m_elec = .200;  
  mass_empty = m_wing + m_htail + m_vtail + m_fuse + m_elec;
end