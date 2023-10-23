%========================================================================
% DESCRIPTION: 
% Matlab script generating figures similar to those in [1].
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
%% Set beam thickness (for normalization)
h_mm = 1;
%% Select Figure from reference [1]
SCENARIO = 'Fig6'; % ['Fig6'|'Fig7'|'Fig8'|'Fig9'|'Fig12a'|'Fig12b']
switch SCENARIO
    case 'Fig6'
        excitationFrequency_Hz = 124;
        baseAcceleration_m_s2 = 14;
        initialSliderPosition_L = 0.3;
        tEnd_s = 60;
        w_h_limit = 2;
    case 'Fig7'
        excitationFrequency_Hz = 104;
        baseAcceleration_m_s2 = 14;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
        w_h_limit = 1;
    case 'Fig8'
        excitationFrequency_Hz = 126;
        baseAcceleration_m_s2 = 10;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
        w_h_limit = 1.5;
    case 'Fig9'
        excitationFrequency_Hz = 124;
        baseAcceleration_m_s2 = 6;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
        w_h_limit = .3;
    case 'Fig12a'
        excitationFrequency_Hz = 106;
        baseAcceleration_m_s2 = 10;
        initialSliderPosition_L = 0.3;
        tEnd_s = 120;
        w_h_limit = 1.5;
    case 'Fig12b'
        excitationFrequency_Hz = 106;
        baseAcceleration_m_s2 = 10;
        initialSliderPosition_L = 0.4;
        tEnd_s = 120;
        w_h_limit = 1.5;
    otherwise
        error(['Unknown scenario ' SCENARIO '.']);
end
%% Load data from files

% sSIM from pseudo-constrained slider (PCS) model
load(['DAT/sim_sSIM_fex' num2str(excitationFrequency_Hz) ...
    '_ExcLevel' num2str(baseAcceleration_m_s2) ...
    '.mat'],'ssSIM_PCS_up','ssSIM_PCS_down',...
    'whatsSIM_PCS_up','whatsSIM_PCS_down');
% Check if there is a jump and truncate if applicable
indJump = find(diff(whatsSIM_PCS_up/(h_mm*1e-3))>.1);
if ~isempty(indJump)
    ssSIM_PCS_up = ssSIM_PCS_up(1:indJump(1));
    whatsSIM_PCS_up = whatsSIM_PCS_up(1:indJump(1));
end
indJump = find(diff(whatsSIM_PCS_down/(h_mm*1e-3))>.1);
if ~isempty(indJump)
    ssSIM_PCS_down = ssSIM_PCS_down(1:indJump(1));
    whatsSIM_PCS_down = whatsSIM_PCS_down(1:indJump(1));
end

% T-w-s data from free slider (FS) model
load(['DAT/sim_FreeSlider_fex' num2str(excitationFrequency_Hz) ...
    '_ExcLevel' num2str(baseAcceleration_m_s2) ...
    '_RelSliderPosStart' strrep(num2str(initialSliderPosition_L),'.','') ...
    '.mat'],'T','w','s');

% T-velocity-s data from experiment
tmp = readmatrix(['DAT/Time_VelBeam_RelSliderPos_' ...
    '_fex' num2str(excitationFrequency_Hz) ...
    '_ExcLevel' num2str(baseAcceleration_m_s2) ...
    '_RelSliderPosStart' strrep(num2str(initialSliderPosition_L),'.','') ...
    '.txt']);
Texp = tmp(:,1);
beamVelocity_m_s = tmp(:,2);
sexp = tmp(:,3);
dtexp = 1e-4; % nominal sampling time in the experiment
%% Extract super-Slow Invariant Manifold (sSIM) from free slider model and 
% experimental data

% Integrate measured velocity to displacement
wexp = cumsum(beamVelocity_m_s*dtexp);
wexp = wexp-smooth(wexp,.1/dtexp);

% Obtain instantaneous amplitude from Hilbert trafo
what = abs(hilbert(w)); 
whatexp = abs(hilbert(wexp)); 

% The super-slow time scale is extracted via moving average
twdw_s = 1;    % window size in s ~ 25 excitation periods
Iwdw = round(twdw_s/(T(2)-T(1))); % windows size in samples
Iwdwexp = round(twdw_s/(Texp(2)-Texp(1))); % windows size in samples

% Take moving average
whatsSIM = smooth(what,Iwdw); 
ssSIM = smooth(s,Iwdw);
whatsSIMexp = smooth(whatexp,Iwdwexp); 
ssSIMexp = smooth(sexp,Iwdwexp);

% Clip off start and end to avoid edge effects
whatsSIM = whatsSIM(Iwdw:end-Iwdw); 
ssSIM = ssSIM(Iwdw:end-Iwdw);
whatsSIMexp = whatsSIMexp(Iwdwexp:end-Iwdwexp); 
ssSIMexp = ssSIMexp(Iwdwexp:end-Iwdwexp);
%% Plot results
figure;
subplot(221); hold on;
title('experiment');
plot(Texp,wexp/(h_mm*1e-3),'k-');  
xlim([0 Texp(end)]); ylim(w_h_limit*[-1 1]);
xlabel('time in s'); ylabel('w(4L/7)/h');
subplot(222); hold on;
title('model');
plot(T,w/(h_mm*1e-3),'g-');  
xlim([0 T(end)]); ylim(w_h_limit*[-1 1]);
xlabel('time in s'); ylabel('w(4L/7)/h');
subplot(223); hold on;
plot(Texp,sexp,'k-','linewidth',2); 
plot(T,s,'g-','linewidth',2); 
xlim([0 T(end)]); ylim([0 .5])
xlabel('time in s'); ylabel('s')
legend('experiment','model','LOCATION','NW');
subplot(224), hold on, box on;
plot(ssSIM_PCS_up,whatsSIM_PCS_up/(h_mm*1e-3),...
    'color',[1 1 1]*.8,'linewidth',3);
plot(ssSIMexp,whatsSIMexp/(h_mm*1e-3),'k-');
plot(ssSIM,whatsSIM/(h_mm*1e-3),'g-');
plot(ssSIM_PCS_down,whatsSIM_PCS_down/(h_mm*1e-3),...
    'color',[1 1 1]*.8,'linewidth',3);
plot(ssSIMexp,whatsSIMexp/(h_mm*1e-3),'k-');
plot(ssSIM,whatsSIM/(h_mm*1e-3),'g-');
plot(ssSIM(1),whatsSIM(1)/(h_mm*1e-3),'go');
plot(ssSIM(end),whatsSIM(end)/(h_mm*1e-3),'gx');
plot(ssSIMexp(1),whatsSIMexp(1)/(h_mm*1e-3),'ko');
plot(ssSIMexp(end),whatsSIMexp(end)/(h_mm*1e-3),'kx');
xlim([0 .5]); ylim([0 2]);
xlabel('s'); ylabel('what(4L/7)/h');
legend('PCS model','experiment','FS model','LOCATION','NW');