function [Comodulogram,p_vec,f_vec] = wave_mod_index_calc(sig_in)
%MOD_INDEX_CALC Summary of this function goes here
%   Detailed explanation goes here
% [P,param_struct] = morlet_transform_fast(s,t,f,fc,FWHM_tc,squared,data_type)

%% Define the amplitude- and phase-frequencies
data_length = size(sig_in,2);
srate = 500; 

phase_fvec = 2:2:50;
amp_fvec = 10:5:125;

phase_fbw = 4;
amp_fbw = 20;


%% Define phase bins
nbin = 18; % number of phase bins
ph_pos=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    ph_pos(j) = -pi+(j-1)*winsize; 
end

%% Filtering and Hilbert transform

Comodulogram = single(zeros(length(phase_fvec),length(amp_fvec)));
trans_ampf = zeros(length(amp_fvec), data_length);
trans_phasef = zeros(length(phase_fvec), data_length);

for ii=1:length(amp_fvec)
    alof = amp_fvec(ii);
    ahif = alof+amp_fbw;
    [~,ampf] = evalc('eegfilt(sig_in,srate,alof,ahif)'); % filtering
    trans_ampf(ii,:) = abs(hilbert(ampf)); % getting the amplitude envelope
end

for jj=1:length(phase_fvec)
    plof = phase_fvec(jj);
    phif = plof + phase_fbw;
    [~,phasef] = evalc('eegfilt(sig_in,srate,plof,phif)'); % filtering 
    trans_phasef(jj, :) = angle(hilbert(phasef)); % getting the phase time series
end


%% Compute MI and comodulogram

cnt1=0;
for ii=1:length(phase_fvec)
    cnt1=cnt1+1;
    %--
    % plof = phase_fvec(ii);
    % phif = plof+phase_fbw;
    %--    
    cnt2=0;
    for jj=1:length(amp_fvec)
        cnt2=cnt2+1;
        %--
        % alof = amp_fvec(jj);
        % ahif = alof+amp_fbw;
        [MI,amp_mu] = mod_index_func_v2(trans_phasef(ii, :), trans_ampf(jj, :), ph_pos);
        Comodulogram(cnt1,cnt2) = MI;
    end
end

p_vec = phase_fvec+phase_fbw/2;
f_vec = amp_fvec+amp_fbw/2;

%% Plot comodulogram
% figure;
% contourf(phase_fvec+phase_fbw/2,amp_fvec+amp_fbw/2,Comodulogram',30,'lines','none');
% set(gca,'fontsize',14);
% ylabel('Amplitude Frequency (Hz)');
% xlabel('Phase Frequency (Hz)');
% colorbar;

% %%  Use the routine below to look at specific pairs of frequency ranges:
% 
% plof = 6;
% phif = 12;
% alof = 60;
% ahif = 100;
% 
% [MI,MeanAmp] = mod_index_fun_v1(lfp,srate,plof,phif,alof,ahif,position);
% 
% bar(10:20:720,[MeanAmp,MeanAmp]/sum(MeanAmp),'k')
% xlim([0 720])
% set(gca,'xtick',0:360:720)
% xlabel('Phase (Deg)')
% ylabel('Amplitude')
% title(['MI = ' num2str(MI)])

end

%% SUBFUNCTIONS
function [mod_index,amp_mu]=mod_index_func_v2(phase_sig, amp_sig, posa)
    % [MI,MeanAmp]=ModIndex_v2(Phase, Amp, position)
    %
    % Phase-amplitude cross-frequency coupling measure:
    %
    % Inputs:
    % Phase = phase time series
    % Amp = amplitude time series
    % position = phase bins (left boundary)
    %
    % Outputs:
    % MI = modulation index (see Tort et al PNAS 2008, 2009 and J Neurophysiol 2010)
    % MeanAmp = amplitude distribution over phase bins (non-normalized)
     

    nbin = length(posa);  
    winsize = 2*pi/nbin;
     
    % now we compute the mean amplitude in each phase:
    amp_mu=zeros(1,nbin); 
    for j=1:nbin   
        ind = (phase_sig < posa(j)+winsize & phase_sig >=  posa(j));
        amp_mu(j) = mean(amp_sig(ind)); 
    end     
    % the center of each bin (for plotting purposes) is position+winsize/2
    
    % quantifying the amount of amp modulation by means of a
    % normalized entropy index (Tort et al PNAS 2008):    
    mod_index=(log(nbin)-(-sum((amp_mu/sum(amp_mu)).*log((amp_mu/sum(amp_mu))))))/log(nbin);
end

%## ==================================================================== %%
function [MI,amp_mu] = mod_index_fun_v1(sig_in,srate,ph_bnd,amp_bnd,posa)
%## 
% Af1 and Af2 define the frequency range investigated as the "amplitude
% modulated" by the phase frequency (e.g., low gamma would be Af1=30 Af2=55)
% 
% position = phase bins (left boundary)
%
% Outputs:
% MI = modulation index
% MeanAmp = Amplitude distribution per phase bin (non-normalized); to
% normalize, do MeanAmp = MeanAmp/sum(MeanAmp)
% the eegfilt routine employed below is obtained from the EEGLAB toolbox 
% (Delorme and Makeig J Neurosci Methods 2004)

%##
Pf1 = ph_bnd(1); %#ok<*NASGU>
Pf2 = ph_bnd(2);
Af1 = amp_bnd(1);
Af2 = amp_bnd(2);

[~,PhaseFreq] = evalc('eegfilt(sig_in,srate,Pf1,Pf2)'); % this is just filtering 
Phase = angle(hilbert(PhaseFreq)); % this is getting the phase time series
[~,AmpFreq] = evalc('eegfilt(sig_in,srate,Af1,Af2)'); % just filtering
Amp = abs(hilbert(AmpFreq)); % getting the amplitude envelope
 
% Now we search for a Phase-Amp relation between these frequencies by
% caclulating the mean amplitude of the AmpFreq in each phase bin of the
% PhaseFreq
 
% Computing the mean amplitude in each phase:

nbin=length(posa);
winsize = 2*pi/nbin;

amp_mu=zeros(1,nbin); 
for j=1:nbin   
    inds = (Phase <  posa(j)+winsize & Phase >=  posa(j));
    amp_mu(j)=mean(Amp(inds)); 
end
% the center of each bin (for plotting purposes) is position+winsize/2
 
% quantifying the amount of amp modulation by means of a
% normalized entropy index (Tort et al PNAS 2008):

MI=(log(nbin)-(-sum((amp_mu/sum(amp_mu)).*log((amp_mu/sum(amp_mu))))))/log(nbin);

end
