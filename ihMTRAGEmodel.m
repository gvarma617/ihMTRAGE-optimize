function [Msing,Mdoub,Mxy] = ihMTRAGEmodel(satdur,dt,rf_t,no_t,b1,delta,...
    pulses,rest_t,alphas,tr,btt,interleave)
% Matlab R2018a code for ihMTRAGEmodel function to simulate magnetization
% from acquisition relating to the inhomogeneous Magnetization Transfer
% (ihMT) experiment.
%
% Input:
% satdur = duration over which to simulate magnetization as it cycles:
%          starting with RF (MT) pulses, then any RF affiliated with on-
%          resonance acquisition (i.e. alphas), then a recovery period
%          (rest_t), and back again. Period = (rf_t+no_t)*pulses + rest_t.
% dt = interval for numerical integration of differential equations
%      relating to two-pool MT model (with and without dipolar reservoirs).
% rf_t = duration over which RF is applied (i.e. MT pulse width) in seconds.
% no_t = duration between successive RF (MT) pulses in units of seconds.
% b1 = B1 amplitude of RF in units of Gauss.
% delta = (off-)resonance frequency of RF (MT) pulse in units of radians.
% pulses = number of RF (MT) pulses to simulate.
% rest_t = duration of acquisition and recovery period in units of seconds.
% alphas = list of flip angles in ° relating to on-resonance acquisition.
% tr = time between any on-resonance pulses as listed in alphas in seconds.
% btt = brain tissue type: 'WM' for white matter; 'GM' for grey matter.
% interleave = logical value relating to sequential or interleaved ihMT
%              acquisition: 0 - sequential, where single and dual frequency
%              experiments remain separated; 1 - interleaved, where the
%              magnetizations from single and dual off-resonance frequency
%              RF (MT) experiments are swapped after Period defined above.
%
% Output:
% Msing = array of magnetizations relating to single frequency MT
%         experiment from free pool, bound pool with dipolar reservoir,
%         bound pool without dipolar reservoir, and dipolar reservoir, for
%         duration satdur in steps of dt.
% Mdoub = array of magnetizations relating to dual frequency MT experiment.
% Mxy = array of transverse magnetizations calculated from free pool
%       components of Msing and Mdoub based on alphas over duration satdur.
%
% Gopal Varma, 2020

global B1RMS
global M_0A
global M_0B1
global M_0B2
global Delta
if nargin < 12
    interleave = 0;
end
if nargin < 11
    btt = 'WM';
end
if strcmp(btt,'WM')
    qihMTparasWM
else
    qihMTparasGM
end
B1RMS = b1;
Delta = delta;

if dt > rf_t || dt > no_t || rf_t/dt < 2 || no_t/dt < 2
    dt = min([rf_t/2,no_t/2]);
    disp(['Time interval changed to ',num2str(dt),'s'])
end
if mod(no_t,dt) > 0
    no_t = round(no_t/dt)*dt;
    disp(['Zero RF time following RF changed to ',num2str(no_t),'s'])
end
if mod(rf_t,dt) > 0
    rf_t = round(rf_t/dt)*dt;
    disp(['Off-resonance RF pulse time changed to ',num2str(rf_t),'s'])
end
if mod(rest_t,dt) > 0
    rest_t = round(rest_t/dt)*dt;
    disp(['Rest time between cycles changed to ',num2str(rest_t),'s'])
end
if nargin > 8
    if mod(tr,dt) > 0
        tr = round(tr/dt)*dt;
        disp(['Time between acquisition RF pulses changed to ',num2str(tr),'s'])
    end
    if rest_t < length(alphas)*tr
        rest_t = length(alphas)*tr;
        disp(['Rest time between cycles (included acquisition) increased to ',...
            num2str(tr),'s'])
    end
end
tot_t = rf_t + no_t;

cycles = ceil(satdur/(tot_t*pulses+rest_t));
Msing = [M_0A,M_0B1,0,M_0B2];
Mdoub = Msing;
Mxy = zeros(cycles * length(alphas), 2);
h = waitbar(0,'Calculating...');
for c = 1:cycles
    for p = 1:pulses
        [~,Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+1):...
            round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+(no_t/dt)+1),:)] = ...
            ode45(@singMT_noRF,0:dt:no_t,...
            Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+1),:));
        [~,Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+(no_t/dt)+1):...
            round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+1),:)] = ...
            ode45(@singMT,no_t:dt:tot_t,...
            Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+(no_t/dt)+1),:));
        [~,Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+1):...
            round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+(no_t/dt)+1),:)] = ...
            ode45(@singMT_noRF,0:dt:no_t,...
            Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+1),:));
        [~,Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+(no_t/dt)+1):...
            round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+1),:)] = ...
            ode45(@dualMT,no_t:dt:tot_t,...
            Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+(p-1)*(tot_t/dt)+(no_t/dt)+1),:));
    end
    if rest_t > 0
        if nargin > 8
            for rep=1:length(alphas)
                [~,Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+(rep-1)*tr/dt+2):...
                    round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+rep*tr/dt+1),:)] = ...
                    ode45(@singMT_noRF,dt:dt:tr,[cos(pi*alphas(rep)/180),1,1,1].*...
                    Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+(rep-1)*tr/dt+1),:));
                Mxy((c-1)*length(alphas)+rep,1) = sin(pi*alphas(rep)/180)*...
                    Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+(rep-1)*tr/dt+1),1);
                [~,Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+(rep-1)*tr/dt+2):...
                    round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+rep*tr/dt+1),:)] = ...
                    ode45(@singMT_noRF,dt:dt:tr,[cos(pi*alphas(rep)/180),1,1,1].*...
                    Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+(rep-1)*tr/dt+1),:));
                Mxy((c-1)*length(alphas)+rep,2) = sin(pi*alphas(rep)/180)*...
                    Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+(rep-1)*tr/dt+1),1);
            end
            if rest_t-length(alphas)*tr > 0
                [~,Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+rep*tr/dt+1):...
                    round(c*(pulses*tot_t+rest_t)/dt+1),:)] = ...
                    ode45(@singMT_noRF,0:dt:rest_t-length(alphas)*tr,...
                    Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+rep*tr/dt+1),:));
                [~,Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+rep*tr/dt+1):...
                    round(c*(pulses*tot_t+rest_t)/dt+1),:)] = ...
                    ode45(@singMT_noRF,0:dt:rest_t-length(alphas)*tr,...
                    Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+rep*tr/dt+1),:));
            end
        else
            [~,Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+1):...
                round(c*(pulses*tot_t+rest_t)/dt+1),:)] = ...
                ode45(@singMT_noRF,0:dt:rest_t,...
                Msing(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+1),:));
            [~,Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+1):...
                round(c*(pulses*tot_t+rest_t)/dt+1),:)] = ...
                ode45(@singMT_noRF,0:dt:rest_t,...
                Mdoub(round((c-1)*((pulses*tot_t+rest_t)/dt)+p*(tot_t/dt)+1),:));
        end
    end
    if interleave
        tmp = Msing(round(c*(pulses*tot_t+rest_t)/dt+1),:);
        Msing(round(c*(pulses*tot_t+rest_t)/dt+1),:) = ...
            Mdoub(round(c*(pulses*tot_t+rest_t)/dt+1),:);
        Mdoub(round(c*(pulses*tot_t+rest_t)/dt+1),:) = tmp;
    end
    waitbar(c/cycles)
end
close(h)

end