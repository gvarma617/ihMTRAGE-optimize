% Matlab R2018a code to create figures relating to optimization of an
% ihMT acquisition within a MPRAGE acquisition, i.e. ihMTRAGE
% Please contact Gopal Varma (gvarma@bidmc.harvard.edu) to report any bugs.

close all, clear variables
set(groot,'defaultLineLineWidth',2,'defaultFigureColor','w',...
    'defaultAxesFontName','Times','defaultAxesFontSize',20)

% Steady-state signal as function of pulse width
pw = 0.0002:0.0002:0.01;
satdur = 5.0; dt = 0.0001; b1 = 0.15; delta = 2*pi*7e3; pulses = 10;
rest_t = 0.0; alphas = 0; tr = 0.0; btt = ['WM';'GM'];
mz_pw = zeros(length(pw), 2, 2);
for p = 1:length(pw)
    for b = 1:2
        [Msing,Mdoub] = ihMTRAGEmodel(5.0,dt,pw(p),19*pw(p),b1,delta,...
            pulses,rest_t,alphas,tr,btt(b,:));
        mz_pw(p,:,b) = [mean(Msing(45001:50000,1),1),mean(Mdoub(45001:50000,1),1)];
    end
end
figure('Name','Steady-state signal as function of pulse width')
plot(1e3*pw,-2*squeeze(diff(mz_pw,1,2))), grid on
xlabel('Pulse width (for fixed 5% duty cycle) [ms]')
ylabel('Simulated ihMT_z')
legend({'WM','GM'},'Location','northwest')

% Build-up with duration of 5% DC preparation
pw = 0.005; no_t = 0.095; dur = 0:dt:satdur; mz_dur = zeros(length(dur),2,2);
for b = 1:2
    [Msing,Mdoub] = ihMTRAGEmodel(satdur,dt,pw,no_t,b1,delta,pulses,...
        rest_t,alphas,tr,btt(b,:));
    mz_dur(:,1,b) = Msing(:,1); mz_dur(:,2,b) = Mdoub(:,1);
end
figure('Name','Build-up with duration of 5% DC preparation')
plot(dur,-2*squeeze(diff(mz_dur,1,2))), grid on
xlabel('RF irradiation duration (5ms every 100ms) [s]')
ylabel('Simulated ihMT_z')

% IhMTz decay during acquisition
fa = [4,8,12]; satdur = 6.0; rest_t = 1.0; alphas = ones(1,125); tr = 0.004;
t = 0:dt:satdur; mz_fa = zeros(length(t),length(fa),2);
for f = 1:length(fa)
    for b = 1:2
        [Msing,Mdoub] = ihMTRAGEmodel(satdur,dt,pw,no_t,b1,delta,pulses,...
            rest_t,fa(f)*alphas,tr,btt(b,:));
        mz_fa(:,f,b) = 2*(Msing(:,1)-Mdoub(:,1));
    end
end
figure('Name','IhMT_z decay during acquisition'), ls = ['- ';'--';': '];
lc = [0,0.447,0.741;0.85,0.325,0.098];
for b = 1:2
    for f = 1:length(fa)
        hold on, plot(t-5,mz_fa(:,f,b),ls(f,:),'Color',lc(b,:))
    end
end
xlabel('Acquisition time following 1 s preparation [s]'), xlim([0 .5])
ylabel('Simulated ihMT_z'), ylim([0 0.25]), grid on
legend({'WM: FA = 4°','          FA = 8°','          FA = 12°',...
    'GM: FA = 4°','         FA = 8°','         FA = 12°'})

% IhMT_z evolution for 5 TR_MTRAGE periods
satdur = 10.0; alphas = 8*ones(1,100); interleave = [0,1];
t = 0:dt:satdur; mz = zeros(length(t),2,2);
for b = 1:2
    for il = 1:length(interleave)
        [Msing,Mdoub] = ihMTRAGEmodel(satdur,dt,pw,no_t,b1,delta,pulses,...
            rest_t,alphas,tr,btt(b,:),interleave(il));
        mz(:,il,b) = 2*(Msing(:,1)-Mdoub(:,1));
    end
end
figure('Name','IhMT_z evolution for 5 TR_MTRAGE periods')
for b = 1:2
    for il = 1:length(interleave)
        hold on, plot(t,mz(:,il,b),ls(il,:),'Color',lc(b,:))
    end
end
xlabel('Scan duration [s]'), grid on
ylabel('Simulated ihMT_z'), ylim([0 0.25])
legend({'WM: sequential','          interleaved','GM: sequential',...
    '         interleaved'})

% IhMT ratio's dependence on acquisition FA
fa = 1:20; satdur = 6.0; no_t = [0.125,0.045]; b1peak = [0.15,0.00]; 
pulses = [8,1]; rest_t = [1.04,0.08]; alphas = [100,10]; pind = [201,401];
MT = zeros(length(fa),2,2,length(b1peak),2);
for f = 1:length(fa)
    for seq = 1:2
        for b1 = 1:length(b1peak)
            for b = 1:2
                [~,~,Mxy] = ihMTRAGEmodel(satdur,dt,pw,no_t(seq),...
                    b1peak(b1),delta,pulses(seq),rest_t(seq),...
                    fa(f)*ones(1,alphas(seq)),tr,btt(b,:));
                MT(f,:,b1,seq,b) = Mxy(pind(seq),:);
            end
        end
    end
end
figure('Name','IhMT ratio''s dependence on acquisition FA')
for b = 1:2
    for seq = 1:2
        hold on, plot(fa,-2*diff(MT(:,:,1,seq,b),1,2)./mean(MT(:,:,2,seq,b),2),...
            ls(seq,:),'Color',lc(b,:))
    end
end
xlabel('Acquisition FA [°]'), grid on
ylabel('Simulated ihMTR'), ylim([0 0.25])
legend({'WM: prepared','          steady state','GM: prepared',...
    '         steady state'})

% IhMT SNR expected after 5 s of each sequence
t = [2.08, 1.3];
figure('Name','IhMT SNR expected after 5 s of each sequence')
for b = 1:2
    for seq = 1:2
        hold on, plot(fa,-2e2*diff(MT(:,:,1,seq,b),1,2)/sqrt(t(seq)),...
            ls(seq,:),'Color',lc(b,:))
    end
end
xlabel('Acquisition FA [°]'), ylim([0 2.5]), grid on
ylabel({'Simulated ihMT_{z} \times sin(FA)';...
    'per square root acquisition time [s^{-1/2}]'})
legend({'WM: prepared','          steady state','GM: prepared',...
    '         steady state'})
