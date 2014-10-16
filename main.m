clear;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RPM =           4000;       % [RPM]
P1 =            214.7700;   % [kPa]
T1 =            54.7;       % [deg C]
Q_f =           20.8;       % [L/h]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_cyl =       4;          % [cylinders]
bore =          0.083;      % [m]
stroke =        0.0904;     % [m]

comp_ratio =    16.5;       % []
vol_eff =       0.92;       % []
comp_eff =      0.7;        % []
turb_eff =      0.9;        % []

LHV =           44000;      % [kJ/kg]
rho_f =         833;        % [kg/m^3]
Q_f =           Q_f/1000/3600/2;    % [m^3/s per cycle]

R =             8.314/29;   % [kJ/kg-K]
n =             1.33;       % []
k =             1.4;        % []

C0 =            0.9328;     % []
C1 =            0;     % []

V(1) =          0;          % Array to store volumes
P(1) =          0;          % Array to store pressures
count =         1;          % Array index counter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Engine Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert T1 to K
T1 = T1 + 273.15;

% V1 = V_displacement/(1-1/r)
V1 = ((bore/2)^2*pi * stroke) / (1 - 1/comp_ratio); % [m^3]

% V2 = V1/r
V2 = V1/comp_ratio;         % [m^3]

% Use ideal gas law to calc mass
m1 = P1*V1/R/T1 * vol_eff;  % [kg]

% Calculate mass in per cycle
m_f = Q_f/(RPM/60) * rho_f; % [kg]

% Calculate heat in per cycle
Q_in = m_f*LHV;             % [kJ]

% Store V and P; Increment counter
V(count) = V1;
P(count) = P1;
count = count + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate P2 for polytropic process
P2 = P1 * comp_ratio^n;     % [m^3]

% Calculate T2 using ideal gas law
T2 = P2*V2/m1/R;            % [K]

% Integrate from V1 to V2 to find compression work
W_comp = 0;                 % Variable to store compression work

% Loop through V1 to V2 to integrate
for Vi=V1:-1e-8:V2
    % Add P(V)dV to W_comp
    W_comp = W_comp + P1*(V1/Vi)^n * -1e-8; % [kJ]
    
    % Store V and P; Increment counter
    V(count) = Vi;
    P(count) = P1*(V1/Vi)^n;
    count = count + 1;
end

W_comp = (P2*V2-P1*V1)/(1-n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Heat I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get internal energy at T2
U2 = m1*getUFromT(T2);      % [kJ]

% Add half of the heat in
U3 = U2 + 0.5*Q_in;         % [kJ]

% Account for heat loss
U3 = U3 - C0*Q_in*RPM^-0.2;

% Add the incoming mass
m3 = m1 + 0.5*m_f;          % [kg]

% Find specific internal energy
u3 = U3 / m3;               % [kJ/kg]

% Volume doesn't change
V3 = V2;                    % [m^3]

% Get T3 from u3
T3 = getTFromU(u3);         % [K]

% Use ideal gas law to find P3
P3 = m3*R*T3/V3;            % [kPa]

% Store V and P; Increment counter
V(count) = V3;
P(count) = P3;
count = count + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Heat II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find enthalpy at T3
H3 = m3*getHFromT(T3);      % [kJ]

% Add half of the heat in
H4 = H3 + 0.5*Q_in;         % [kJ]

% Account for heat loss
H4 = H4 - C0*Q_in*RPM^-0.2;

% Add the incoming mass
m4 = m3 + 0.5*m_f;          % [kg]

% Find specific enthalpy
h4 = H4 / m4;               % [kJ/kg]

% Pressure doesn't change
P4 = P3;                    % [kPa]

% Get T4 from h4
T4 = getTFromH(h4);         % [K]

% Find V4 using ideal gas law
V4 = m4*R*T4/P4;            % [m^3]

% Calculate PdV work
W_heat = P3 * (V4 - V3);

% Store V and P; Increment counter
V(count) = V4;
P(count) = P4;
count = count + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Final volume = initial volume
V5 = V1;                    % [m^3]

% Calculate P5 using polytropic expansion
P5 = P4 * (V4/V5)^n;        % [kPa]

% Calculate T5 using ideal gas law
T5 = P5*V5/m4/R;            % [K]

% Calculate internal energy at T5
U5 = m4*getUFromT(T5);

% Integrate from V5 to V1 to find expansion work
W_exp = 0;                  % Variable to store expansion work

% Loop through V4 to V5 to integrate
for Vi=V4:1e-8:V5
    W_exp = W_exp + P4*(V4/Vi)^n * 1e-8;    % [kJ]
    
    % Store V and P; Increment counter
    V(count) = Vi;
    P(count) = P4*(V4/Vi)^n;
    count = count + 1;
end

W_exp = (P5*V5-P4*V4)/(1-n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exhaust
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cacculate internal energy at T1
U1 = m4*getUFromT(T1);

% Store V and P; Increment counter
V(count) = V1;
P(count) = P1;
count = count + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot P-V diagram
plot(V, P);

% Output brake MEP
BMEP = (W_exp + W_heat + W_comp) / (V1-V2) / 100

% Output net power calculation
dWdt = (W_exp + W_heat + W_comp) * RPM / 60 * 4     % [kW]

% Calculate friction losses
Sbar_p = RPM/60*2*stroke;
PME = 75 + 48*RPM/1000+0.4*Sbar_p^2;
dWdt_fric = C1*PME*((bore/2)^2*pi * stroke)*RPM/2

% Output brake power calculation
dWdt_brake = dWdt - dWdt_fric