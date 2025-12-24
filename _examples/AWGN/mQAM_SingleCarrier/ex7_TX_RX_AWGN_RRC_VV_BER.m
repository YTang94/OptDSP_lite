% 16QAM + AWGN + RRC pulse shaping/matched filter + VV CPE + BER

clear;
close all;
clc;

% Add OptDSP library to path
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(repoRoot));

%% Simulation Parameters
Fs = 80e9;              % sampling rate [Sa/s]
Rs = 10e9;              % symbol rate [Baud]
nSpS = Fs/Rs;           % samples per symbol (should be integer)
nSyms = 2^14;           % number of symbols
nSamples = nSyms * nSpS;
rollOff = 0.1;
SNR_dB = 20;
SNR_sweep_dB = 10:2:24;

PARAM = setSimulationParams(Fs,nSamples);

SIG.symRate = Rs;
SIG.M = 16;
SIG.nPol = 1;
SIG.rollOff = rollOff;
SIG = setSignalParams(SIG,PARAM);

QAM = QAM_config(SIG);

%% Generate Bits and QAM Symbols
txBits = randi([0 1],SIG.nPol,SIG.nBits) > 0;
[Stx_syms,txSyms] = Tx_QAM(QAM,txBits);

%% Pulse Shaping (RRC)
PS.type = 'RRC';
PS.rollOff = rollOff;
[Stx,PS] = pulseShaper(Stx_syms,nSpS,PS);

%% AWGN Channel
[Srx,~,SNR_out] = setSNR(Stx,SNR_dB,Fs,Rs);
fprintf('Target SNR = %.2f dB, Measured SNR = %.2f dB\n',SNR_dB,SNR_out);

%% Matched Filter (RRC)
LPF.type = 'RRC';
LPF.rollOff = rollOff;
[Srx,LPF] = LPF_apply(Srx,LPF,Fs,Rs);

%% Carrier Phase Estimation (Viterbi & Viterbi)
CPE = presetCPE('method','VV','mQAM','16QAM','nTaps',64,'nSpS',nSpS);
CPE.decision = 'QPSKpartition';
[Srx,CPE] = carrierPhaseEstimation(Srx,Stx_syms,CPE,QAM.IQmap,nSpS);

%% Downsample to 1 Sample/Symbol
Srx_syms = Srx(:,1:nSpS:end);

%% Synchronize Rx and Tx Symbol Streams
SYNC.method = 'complexField';
[Srx_sync,delay] = syncSignals(Stx_syms(1,:),Srx_syms(1,:),SYNC);
fprintf('Estimated symbol delay = %d\n',delay);

%% Symbol Decision and BER
rxSymInd = signal2symbol(Srx_sync,QAM.IQmap,[]);
rxBits = sym2bit(rxSymInd,QAM.nBpS);
[BER,errPos] = BER_eval(txBits,rxBits);
fprintf('BER = %.3e (errors: %d)\n',BER,numel(errPos));

%% Visualizations
% Constellation (Tx and Rx after CPE)
figure('Name','Constellation');
subplot(1,2,1);
plot(real(Stx_syms(1,1:5000)),imag(Stx_syms(1,1:5000)),'.');
axis equal; grid on;
title('Tx Constellation (Symbols)');
xlabel('I'); ylabel('Q');

subplot(1,2,2);
plot(real(Srx_sync(1,1:5000)),imag(Srx_sync(1,1:5000)),'.');
axis equal; grid on;
title('Rx Constellation (After CPE)');
xlabel('I'); ylabel('Q');

% Time-domain waveform (I/Q)
figure('Name','Time-Domain Waveform');
t = (0:size(Stx,2)-1)/Fs;
plot(t(1:2000),real(Stx(1,1:2000)),'-');
hold on;
plot(t(1:2000),imag(Stx(1,1:2000)),'-');
grid on;
title('Tx Waveform (RRC Shaped)');
xlabel('Time [s]');
ylabel('Amplitude');
legend('I','Q');

% BER vs SNR curve
BER_curve = zeros(size(SNR_sweep_dB));
for k = 1:numel(SNR_sweep_dB)
    [Srx_k,~,~] = setSNR(Stx,SNR_sweep_dB(k),Fs,Rs);
    [Srx_k,~] = LPF_apply(Srx_k,LPF,Fs,Rs);
    [Srx_k,~] = carrierPhaseEstimation(Srx_k,Stx_syms,CPE,QAM.IQmap,nSpS);
    Srx_syms_k = Srx_k(:,1:nSpS:end);
    Srx_sync_k = syncSignals(Stx_syms(1,:),Srx_syms_k(1,:),SYNC);
    rxSymInd_k = signal2symbol(Srx_sync_k,QAM.IQmap,[]);
    rxBits_k = sym2bit(rxSymInd_k,QAM.nBpS);
    BER_curve(k) = BER_eval(txBits,rxBits_k);
end

figure('Name','BER Curve');
semilogy(SNR_sweep_dB,BER_curve,'-o');
grid on;
xlabel('SNR [dB]');
ylabel('BER');
title('BER vs SNR (16QAM, RRC, VV CPE)');
