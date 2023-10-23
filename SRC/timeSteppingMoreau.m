%========================================================================
% DESCRIPTION: 
% Function carrying out time step integration using Moreau's scheme. 
% The considered measure differential inclusion is
%       M*du + h(q,u,t)*dt + W*dP = 0.
% Contact activity and kinematics is defined in the given function. 
% At each contact, unilateral and dry frictional interactions are 
% modeled by the Signorini and Coulomb laws combined with Newton’s impact 
% law. The dynamic contact problem is resolved using an Augmented
% Lagrangian approach.
% 
% INPUT
%   ts      start time
%   te      end time
%   q0      initial coordinates
%   u0      initial velocities
%   M       mass matrix (assumed constant in this code)
%   h       function handle for remaining forces
%   conKin  function handle for contact activity and kinematics 
%               (t,q) --> (IC, W, wt, gN) stored in structure con.
%   nCmax   total number of potential contacts
%   eN      restitution coefficient in normal direction (eT=0)
%   mu      friction coefficient
%   dtTMP   target time step size
%   rtol    relative tolerance (convergence of contact percussions)
%   atol    absolute tolerance (...)
% 
% OUTPUT
%   T       row vector of time levels
%   Q       coordinates, dim(q0) x dim(T)
%   U       velocities, dim(u0) x dim(T)
%   PP      contact percussions
%   GN      normal gaps
%   ACT     contact activity
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
function [T,Q,U,PP,GN,ACT] = timeSteppingMoreau(ts,te,q0,u0,...
    M,h,conKin,nCmax,eN,mu,...
    dtTMP,rtol,atol,displayProgress)
%% Handle input, define contact settings and initialize output

% Reset time step size to have equal steps from start to end
Ntd = round((te-ts)/dtTMP);
dt = (te-ts)/Ntd;
T = ts:dt:te;

% Contact settings
nCdim = 2;
projC = @(xi) proxUF(xi(1),xi(2),mu);
epsI = [eN;0]; % second entry eT=0

% Initialize output data
n = length(q0);
Q = zeros(n,Ntd+1);
U = zeros(n,Ntd+1);
% additional output on request: 
if nargout>3
    ACT = zeros(nCmax,Ntd+1);       % contact activity
    PP = zeros(nCdim*nCmax,Ntd+1);  % contact percussions
    GN = NaN*ones(nCmax,Ntd+1);     % normal gaps (known only for active)
end

% Store initial values
Q(:,1) = q0; 
U(:,1) = u0;
% It is assumed that no contact is active at time ts!
P(1:nCmax) = {zeros(nCdim,1)};
%% Loop over time levels
for j=1:Ntd
    % Time and states at beginning
    tB = T(j); 
    qB = Q(:,j); 
    uB = U(:,j);
    
    % Values at mid-point (Euler half step)
    tM = tB+dt/2; 
    qM = qB+dt/2*uB;
    hM = h(tM,qM,uB);
    
    % Evaluate contact kinematics at mid-point
    con = conKin(tM,qM);
    
    % Determine active set and corresponding kinematics
    IC = con.IC;
    nC = length(IC);
    W = con.W; 
    wt = con.wt;
    
    % If at least one contact is active, solve contact problem
    if nC>0
        %% Setup inclusion problem
        
        % Determine contact force directions
        WM = zeros(n,nCdim*nC); wM = zeros(nCdim*nC,1);
        for ii=1:length(IC)
            i = IC(ii);
            II = nCdim*(ii-1)+(1:nCdim);
            WM(:,II) = W{i}(tM,qM);
            wM(II) = wt{i}(tM,qM);
        end
        
        % Determine Delassus matrix G and vector c
        G = WM'*(M\WM);
        c = kron( eye(nC), eye(nCdim) + diag(epsI) )*WM'*uB + ...
            WM'*(M\(hM*dt)) + ...
            kron( eye(nC), eye(nCdim) + diag(epsI) )*wM;
        
        % Define Augmented Lagrangian parameter
        r_al = zeros(nC,1);
        G_strictly_diagonal_dominant = ...
            all((2*abs(diag(G))) >= sum(abs(G),2));
        for ii=1:nC
            II = nCdim*(ii-1)+(1:nCdim);
            if G_strictly_diagonal_dominant
                r_al(ii) = 1/max(diag(G(II,II)));
            else
                r_al(ii) = 1/max(sum(abs(G(II,:)),2));
            end
        end
        %% Iterative solution of the inclusion problem
        converged = 0;          % Boolean flag indicating convergence
        Pk = P(IC);             % adopt previous percussions
        Pkpi = cell(size(Pk));  % initialize percussion cell array
        k = 0;                  % iteration counter
        while ~converged
            % Apply successive over-relaxation
            err = 0;
            for ii=1:nC
                II = nCdim*(ii-1)+(1:nCdim);
                Pkpi{ii} = projC( Pk{ii} - r_al(ii)*...
                    (G(II,:)*vertcat(Pk{:})+c(II)) );
                err = err + norm(Pkpi{ii}-Pk{ii});
                Pk{ii} = Pkpi{ii};
            end
            % Evaluate convergence criterion
            converged = err<=norm(vertcat(Pk{:}))*rtol+atol;
            k = k + 1;
        end
    end
    
    % Determine states at the end of the time step tE
    tE = tM + dt/2;
    if nC==0
        uE = uB + M\( hM*dt );
    else
        uE = uB + M\( hM*dt + WM*vertcat(Pk{:}) );
    end
    qE = qM + dt/2*uE;
    
    % Store output quantities
    T(j+1) = tE; 
    Q(:,j+1) = qE; 
    U(:,j+1) = uE;
    
    if nC>0
        P(IC) = Pk;
        % additional output on request
        if nargout>3
            GN(:,j+1) = con.gN;
            ACT(con.IC,j+1) =  con.IC;
            ICC = kron(nCdim*(IC-1),ones(nCdim,1)) + ...
                repmat((1:nCdim)',nC,1);
            PP(ICC,j+1) = vertcat(Pk{:});
        end
    end

    % Report status every percent of simulation time
    if (nargin<=13 || displayProgress) && ...
            (mod(j,round(Ntd/100))==0 || j==Ntd)
        disp(['Simulation at ' num2str(round(j/Ntd*100)) '%.']);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Proximal point function for 2D unilateral-frictional contact law
function prxi = proxUF(xiN,xiT,mu)
    prxi = zeros(2,1);
    % Unilateral contact law
    prxi(1) = max(xiN,0);
    % Frictional contact law
    if xiT>mu*prxi(1)
        prxi(2) = mu*prxi(1);
    elseif xiT<-mu*prxi(1)
        prxi(2) = -mu*prxi(1);
    else
        prxi(2) = xiT;
    end
end