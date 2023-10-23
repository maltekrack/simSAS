%========================================================================
% DESCRIPTION: 
% This function identifies the active set of contacts and required 
% kinematic quantities for the self-adaptive beam-slider system.
% Some simplifications are already build into the code, incl. the small
% curvature approximation dw^2<<1, assumption of super-slow time scale of
% transport along beam. Those simplifications were found to be valid, based
% on the excellent agreement with the results obtained for the initial
% model, where those simplifications were not made.
% 
% INPUT
%   t       time
%   q       vector of generalized coordinates with q = [qb;qs] and 
%           qs = [xC;zC;beta]. C stands for center of mass.
%   PHI     function handle for the beam's modal deflection shape in z dir.
%   dPHI    function handle for the beam's modal slpe shape in z direction
%   w0      function handle for base motion in z direction
%   dw0     function handle for base velocity in z direction
%   L       length of beam
%   tb      thickness of beam
%   rCP_S   vector from slider's center of mass to contact points in slider
%           frame of reference, dimension 2 (plane) x 4 (4 contact points)
%   zetaZ   vertical indicator for contact points
%   rtol    tolerance for contact activation, relative to tb2
% 
% OUTPUT
%   con.IC  contact activity index
%   con.gN  normal gap (>0 means inactive)
%   con.W   contact force directions (dg/dq)'
%   con.wt  imposed component of gap velocity dg/dt
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
function con = contactBeamSlider(t,q,PHI,dPHI,w0,dw0,L,tb,rCP_S,zetaZ,...
    varargin)

% Set tolerance for contact activation, rel. to beam half thickness
rtol = 1e-5;    % default value
if nargin>10 && ~isempty(varargin{11})
    rtol = varargin{11}; % use specified value
end
tb2 = tb/2; % beam half thickness

% Determine beam 

% Initialize output
nC = size(rCP_S,2); % number of potential contacts (=4 in our case)
IC = zeros(nC,1);
W = cell(nC,1);
wt = W;
gN = zeros(size(W));

% Extract beam and slider coordinates
qb = q(1:end-3);
qs = q(end-2:end);
xC = qs(1); zC = qs(2); beta = qs(3);

% Rotation matrix to inertial (N) from slider (S) frame of reference
A_NS = [cos(beta) sin(beta); -sin(beta) cos(beta)];

% Position of slider's center of mass (C) and its contact points (P)
rC_N = [xC;zC];
rP_N = repmat(rC_N,1,nC) + A_NS*rCP_S;

% Loop over potential contacts
for ip=1:nC
    xPi = rP_N(1,ip);
    zPi = rP_N(2,ip);
    zB_i = w0(t) + PHI(xPi/L)*qb;
    if zetaZ(ip)*(zPi-zB_i)-tb2<tb2*rtol
        %% Contact possible, check for actual contact
        
        % Normalized contact location
        si = xPi/L;
        
        % Beam's elastic transversal displacement and slope
        wi = PHI(si)*qb;
        dwi = dPHI(si)*qb;
        
        % Derivative of slider rotation matrix and contact point location
        dA_NS_dbeta = [0 1;-1 0]*A_NS;
        drP_dqs = [eye(2) dA_NS_dbeta*rCP_S(:,ip)];
        
        % Rotation angle of beam, rotation matrix to contact frame of ref.
        alp = atan(dwi);
        dalp_dqb = dPHI(si);
        A_CN = zetaZ(ip)*[-sin(alp) cos(alp); -cos(alp) -sin(alp)];
        dA_CN_dalp= [0 1;-1 0]*A_CN;
        
        % Location of projected point on beam, and derivatives
        rB_proj = [xPi;w0(t)+wi] + zetaZ(ip)*tb2*[-dwi;1];
        drB_proj_dqb = [zeros(1,size(qb,1)); PHI(si)] + ...
            zetaZ(ip)*tb2*( [-dPHI(si);zeros(1,size(qb,1))] - ...
            [-dwi^2;dwi]*dPHI(si) );
        drB_proj_dt = [0;dw0(t)];
        
        % Gap, and derivatives
        g = A_CN*( rP_N(:,ip) - rB_proj );
        dg_dqs = A_CN*( drP_dqs );
        dg_dqb = -A_CN*drB_proj_dqb + ...
            dA_CN_dalp*( rP_N(:,ip) - rB_proj )*dalp_dqb;
        dg_dq = [dg_dqb dg_dqs];
        dg_dt = @(t) -A_CN*( drB_proj_dt );
        
        % Determine sign-valued gap and contact activity
        gN(ip) = g(1);
        IC(ip) = gN(ip)<0;
        
        % Define contact force directions and imposed gap velocity
        W{ip} = @(t,q) dg_dq';
        wt{ip} = @(t,q) dg_dt(t)';
    end
end

% Colled output in structure
con.IC = find(IC);
con.W = W;
con.wt = wt;
con.gN = gN;