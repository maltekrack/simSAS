%========================================================================
% DESCRIPTION: 
% Matlab script simulating the behavior of a certain self-adaptive system (SAS).
%========================================================================
% This file is part of simSAS.
% 
% If you use simSAS, please refer to either one or mulitple of these
% articles:
%   [1] Müller, F.; Beck, M.; Krack, M.: Experimental validation of a model 
%       for a self-adaptive beam–slider system. MSSP, 2023. 
%       doi: 10.1016/j.ymssp.2022.109551
%       MODEL EXTENDED BY IMPERFECT BOUNDARY CONDITIONS &
%       EXPERIMENTAL VALIDATION
%   [2] Müller, F.; Krack, M.: Explanation of the self-adaptive dynamics of 
%       a harmonically forced beam with a sliding mass. AAM, 2020. 
%       doi: 10.1007/s00419-020-01684-5
%       MODEL EXTENDED BY GEOMETRIC NONLINEARITY OF THE BEAM
%   [3] Krack, M.; et al.: Toward understanding the self-adaptive dynamics 
%       of a harmonically forced beam with a sliding mass. AAM, 2017. 
%       doi: 10.1007/s00419-016-1218-5
%       INITIAL MODEL
% 
% COPYRIGHT AND LICENSING: 
% Copyright (C) 2023  Malte Krack (malte.krack@ila.uni-stuttgart.de) and
%                     Florian Müller (florianm987@gmail.com)
% This program comes with ABSOLUTELY NO WARRANTY. 
% simSAS is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
clearvars;
close all;
clc;
addpath('SRC');
%% System parameters 
L_mm = 140;         % free length (beam)
h_mm = 1;           % thickness (beam)
rho_kg_m3 = 7683;   % density (beam)
rhoAL_kg = 15.1e-3; % mass (beam, free length)
E_Pa = 210e9;       % Young's modulus (beam)
m_kg = 46.2e-3;     % mass (slider)
JC_kgmm2 = 3.15;    % rotary inertia w.r.t. center of mass (slider)
B_mm = 10;          % distance between contacts (slider, on same side)
d_mm = 4.1;         % distance geometric center to center of mass (slider)
R_h = 1.05;         % relative gap size
kt_N_m = 1.93e7;    % axial clamping stiffness
kr_Nm_rad = 125;    % rotational clamping stiffness
Dmod = 0.1e-2;      % modal damping ratio (beam; set equal for all modes)
mu = 0.2;           % friction coefficient (beam-slider)
eN = .5;            % restitution coefficient in normal direction (eT=0)
%% Excitation and initial conditions, simulation end time

% Select scenario, labeled after Figures in reference [1]
SCENARIO = 'Fig6'; % ['Fig6'|'Fig7'|'Fig8'|'Fig9'|'Fig12a'|'Fig12b']

switch SCENARIO
    case 'Fig6'
        % SIGNATURE MOVE
        excitationFrequency_Hz = 124;
        baseAcceleration_m_s2 = 14;
        initialSliderPosition_L = 0.3;
        tEnd_s = 60;
    case 'Fig7'
        % TRIVIAL ADAPTION
        excitationFrequency_Hz = 104;
        baseAcceleration_m_s2 = 14;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
    case 'Fig8'
        % RELAXATION OSCILLATION
        excitationFrequency_Hz = 126;
        baseAcceleration_m_s2 = 10;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
    case 'Fig9'
        % SELF-DAMPING
        excitationFrequency_Hz = 124;
        baseAcceleration_m_s2 = 6;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
    case 'Fig12a'
        % TRIVIAL ADAPTION
        excitationFrequency_Hz = 106;
        baseAcceleration_m_s2 = 10;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
     case 'Fig12b'
        % SIGNATURE MOVE
        excitationFrequency_Hz = 106;
        baseAcceleration_m_s2 = 10;
        initialSliderPosition_L = 0.4;
        tEnd_s = 120;
    otherwise
        error(['Unknown scenario ' SCENARIO '.']);
end
%% Slider geometry: vectors from center of mass to potential contacts

% Inner clearance
R_mm= R_h*h_mm;

% Location of center of mass (C)
rC_S = [0;-d_mm];

% Vertical (zetaZ) and horizontal (zetaX) indices of slider contact points
% NOTE: We count counter-clockwise from top-left.
zetaZ = [1 -1 -1 1]; % 1: top, -1: bottom
zetaX = [-1 -1 1 1]; % 1: right, -1: left

% Vectors from the center of mass to the four slider contact points
rCP_S = zeros(size(rC_S,1),size(zetaZ,2));
for ip=1:size(rCP_S,2)
    rPip_S = [zetaX(ip)*B_mm/2;zetaZ(ip)*R_mm/2];
    rCP_S(:,ip) = rPip_S - rC_S;
end

% Convert from mm to m
rC_S = rC_S*1e-3;
rCP_S = rCP_S*1e-3;
%% Set up modal beam model

% Auxiliary parameters of the cross section
width_mm = rhoAL_kg/(L_mm*1e-3)/rho_kg_m3/(h_mm*1e-3)*1e3;
I_m4 = 1/12*width_mm*h_mm^3*1e-12;
A_m2 = width_mm*h_mm*1e-6;

% Modal truncation order
Nmod = 5;

% Call function that sets up the actual model, see function for definition
% of output.
[Mb,hb,PHI,dPHI] = modalBeamModel(L_mm*1e-3,rho_kg_m3,E_Pa,...
    A_m2,I_m4,kr_Nm_rad,kt_N_m,Nmod,Dmod);
%% Set up system model (beam and slider, base motion)
% The measure differential inclusion takes the form:
%       M*du + h(q,u,t)*dt + W*dP = 0.
% NOTE: We work in SI units.

% Harmonic base motion
Om = 2*pi*excitationFrequency_Hz;
w0c = baseAcceleration_m_s2/Om^2;
w0 = @(t) w0c*cos(Om*t);
dw0 = @(t) -w0c*Om*sin(Om*t);
ddw0 = @(t) -w0c*Om^2*cos(Om*t);

% Mass matrix
Ms = diag([m_kg m_kg JC_kgmm2*1e-6]);   % slider mass matrix
M = blkdiag(Mb,Ms);

% Internal and external forces (except inertia and contact)
gacc_m_s2 = 9.81;                       % gravitational acceleration
hs = [0;-m_kg*gacc_m_s2;0];             % slider forces
h = @(t,q,u) [hb(q(1:Nmod),u(1:Nmod),ddw0(t));hs];

% Contact activity and kinematics
contactActivityAndKinematics = @(t,q) contactBeamSlider(t,q,PHI,dPHI,...
    w0,dw0,L_mm*1e-3,h_mm*1e-3,rCP_S,zetaZ);
nCmax = 4; % total number of potential contacts
%% Simulation

% Set simulation start and end times
tStart = 0;
tEnd = tEnd_s;

% Initial values: slider at given xC, contacts sym. open, otherwise zero
q0 = zeros(Nmod+length(Ms),1);
q0(end-2) = initialSliderPosition_L*L_mm*1e-3;
q0(end-1) = rC_S(2)+w0(tStart);
u0 = zeros(Nmod+length(Ms),1);

% Time step size and solver tolerances
dt = 2e-5; % ~2*pi/om(5)/14
rtol = 1e-5; atol = 1e-7; % relative and absolute tol. for impl. proj.

% Time step integration using Moreau-type scheme
[TT,QQ,UU] = timeSteppingMoreau(tStart,tEnd,q0,u0,...
    M,h,contactActivityAndKinematics,nCmax,eN,mu,...
    dt,rtol,atol);
%% Reduce and store data

% Take only every dNpt point
dNpt = 10;
T = TT(1:dNpt:end); 
w = PHI(4/7)*QQ(1:Nmod,1:dNpt:end);   % elastic beam displacement at 4/7*L
s = QQ(end-2,1:dNpt:end)/(L_mm*1e-3); % relative slider location xC/L

% Save data to file
save(['DAT/sim_FreeSlider_fex' num2str(excitationFrequency_Hz) ...
    '_ExcLevel' num2str(baseAcceleration_m_s2) ...
    '_RelSliderPosStart' strrep(num2str(initialSliderPosition_L),'.','') ...
    '.mat'],'T','w','s');
