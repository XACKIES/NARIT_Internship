%% ADS-B Receiver PlutoSDR ( Model 2 : Frame Buffer )
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s) 
%     -> FIR pre-filter 
%     -> magnitude-squared (energy)
%     -> Frame Buffer ( 10 Frame of ADS-B Message : 1.2ms/ Sliding buffer length )
%     -> matched filter (preamble template)
%     -> correlation peak search 
%     -> preamble validation (chip-level energy check)
%     -> PPM demodulation (early-late method)
%     -> CRC-24 parity check
%
% References:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS).
%   [2] ICAO Annex 10, Vol. IV.
%   [3] CRC polynomial: 0x1FFF409.
%   [4] Oppenheim & Schafer, DSP (Linear-phase FIR).
 
clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT   = false;      % Enable/disable debug plots
DEBUG_EVERY  = 10;         % Plot every N frames
iterCount    = 0;

%% ---------- Radio / Buffer ----------
fc       = 1090e6;         % ADS-B carrier frequency
sampRate = 10e6;           % 10 MS/s -> 5 samples per 0.5 us chip at 2 MHz chip rate
frameLen = 65536;          % SDR read size (samples per hardware read)
bufferLen= 16*frameLen;    % Sliding buffer length (for correlation/search windowing)

% PlutoSDR front-end with Fast Attack AGC
rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...f  +
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...   % AGC enabled for robust input level
  'OutputDataType','double');

disp('✅ PlutoSDR ready (AGC Fast Attack) ... Listening 1090 MHz ...');

%% ---------- FIR front-end (pre-buffer) ----------
% User-provided linear-phase FIR (real coeffs), applied to complex IQ.
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286 -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];

% Streaming filter state (kept across frames)
firState = zeros(length(b_fir)-1, 1);

%% ---------- ADS-B PHY Parameter ----------
spc_chip  = round(sampRate/2e6);  % Samples per chip (chip = 0.5 us)
spb_bit   = 2*spc_chip;           % Samples per bit (bit = 1.0 us)

% 16-chip preamble (8 us) in chip-domain (1 = pulse present, 0 = gap)
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];

ADS_B_Parameter.SamplesPerChip             = spc_chip;
ADS_B_Parameter.SyncSequence               = SyncSequence;
ADS_B_Parameter.SyncSequenceLength         = length(SyncSequence);
ADS_B_Parameter.SyncSequenceHighIndices    = find(SyncSequence==1);
ADS_B_Parameter.SyncSequenceLowIndices     = find(SyncSequence==0);
ADS_B_Parameter.SyncSequenceNumHighValues  = numel(ADS_B_Parameter.SyncSequenceHighIndices);
ADS_B_Parameter.SyncSequenceNumLowValues   = numel(ADS_B_Parameter.SyncSequenceLowIndices);
ADS_B_Parameter.PreambleLength             = ADS_B_Parameter.SyncSequenceLength * spc_chip;   % 16 * spc
ADS_B_Parameter.LongPacketLength           = 112 * spb_bit;                             % 112 bits * samples/bit
ADS_B_Parameter.MaxPacketLength            = ADS_B_Parameter.PreambleLength + ADS_B_Parameter.LongPacketLength;
ADS_B_Parameter.MaxNumPacketsInFrame       = 64;
ADS_B_Parameter.SyncDownsampleFactor       = 1;

conting = 0;
% Matched-filter template : ±1 version of the preamble replicated at chip rate
preamble_bip = 2*SyncSequence - 1;                % Map {0,1} -> {-1,+1}
mf = flipud(repelem(preamble_bip(:), spc_chip));  % Time-reversed template (length = 16*spc)

% Sliding buffer for streaming correlation
xBuff = zeros(bufferLen,1);

%% ---------- Main loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Read RF samples
  newFrame = rx();
  if isempty(newFrame), continue; end

  % 2) FIR pre-filter (streaming with state); apply to complex IQ
  [newFrameFIR, firState] = filter(b_fir, 1, newFrame, firState);

  % 3) Update sliding buffer with FIR'ed samples
  xBuff = [xBuff(numel(newFrameFIR)+1:end); newFrameFIR];

  % 4) Convert to instantaneous energy |u|^2
  energySig = abs(xBuff).^2;

  % 5) Matched-filter correlation against the preamble (on the energy signal)
  xFilt = conv(energySig, mf, 'same');

  % 6) Framework-like packet search on the correlation output
  [packetSamples, packetCnt, syncTimeVec] = packetSearch(abs(xFilt), xBuff, energySig, ADS_B_Parameter);

  % 7) Optional debug plots (every DEBUG_EVERY frames)
  if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY)==0
    try
      figure(1); clf;
      t = (1:numel(energySig))/sampRate*1e3; % Time axis in ms

      subplot(2,1,1);
      plot(t, energySig); grid on;
      title('Energy |u|^2 (buffer)'); xlabel('Time (ms)'); ylabel('Energy');
      hold on;
      if ~isempty(syncTimeVec)
        plot(t(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), ...
             energySig(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), 'ro');
        legend('Energy','Sync candidates');
      end

      subplot(2,1,2);
      plot(t, xFilt); grid on;
      title('Matched-filter correlation'); xlabel('Time (ms)'); ylabel('Correlation');
      hold on;
      if ~isempty(syncTimeVec)
        plot(t(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), ...
             xFilt(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), 'rx');
        legend('Correlation','Sync candidates');
      end
      drawnow limitrate;
    catch
      end
  end

  % 8) Packet decode
  for k = 1:packetCnt
    region = abs(packetSamples(:,k));

    % Safety: reshape region into [samples-per-bit x num-bits]
    nBits = floor(numel(region)/spb_bit);
    if nBits < 112, continue; end
    segs  = reshape(region(1:nBits*spb_bit), spb_bit, nBits);

    % PPM demodulation: compare early (0–0.5 us) vs late (0.5–1.0 us) energy
    early = sum(segs(1:spc_chip, :), 1);
    late  = sum(segs(spc_chip+1:end, :), 1);
    bits  = double(early > late);
    bits  = bits(1:112);   % Keep the standard 112-bit payload (DF17)

    % CRC-24 parity check across first 112 bits (DO-260 Mode S parity)
    if checkCRC(bits)
      DF = bin2dec(num2str(bits(1:5)));
      if DF==17
        conting = conting+1;
        CA   = bin2dec(num2str(bits(6:8)));
        ICAO = bits(9:32);
        DATA = bits(33:88);
        disp('========= VALID ADS-B (DF=17) =========');
       
      fprintf('CA   = %d\n', CA);
      fprintf('ICAO = %06X\n', bin2dec(char(ICAO + '0')));
      fprintf('DATA = %s\n', char(DATA + '0'));
      fprintf('#Aircraft Detecting Counting = %d\n', (conting));
      disp('======================================');
      
      end
    end
  end
end

%% -------- helper: Framework-like packet search --------
function [packetSamples, packetCnt, syncTimeVec] = packetSearch(xFilt, xBuff, energySig, ADS_B_Parameter)
  % Inputs:
  %   xFilt     - correlation output (abs / correlation magnitude)
  %   xBuff     - current sliding buffer of raw IQ samples
  %   energySig - |u|^2 energy sequence
  %   ADS_B_Parameter - struct with PHY parameters 

  % Outputs:
  %   packetSamples - [LongPacketLength x MaxNumPacketsInFrame] extracted regions (data only)
  %   packetCnt     - number of packets found in this frame
  %   syncTimeVec   - sync start indices (buffer positions) for found packets

  spc = ADS_B_Parameter.SamplesPerChip;
  syncLen    = ADS_B_Parameter.SyncSequenceLength;
  syncSigLen = syncLen*spc;
  xLen       = length(xBuff);

  % Subframe length ~ maximum packet length (preamble + 112 bits)
  subFrameLen     = ADS_B_Parameter.MaxPacketLength;
  subFrameDownLen = subFrameLen / ADS_B_Parameter.SyncDownsampleFactor;
  numSubFrames    = floor(xLen / subFrameLen);

  packetSamples = zeros(ADS_B_Parameter.LongPacketLength, ADS_B_Parameter.MaxNumPacketsInFrame, 'like', xBuff);
  syncTimeVec   = zeros(ADS_B_Parameter.MaxNumPacketsInFrame,1);
  packetCnt     = 0;

  for p = 0:(numSubFrames-2)
    % (A) Find correlation peak within this subframe
    idx = double(p)*subFrameDownLen + (1:subFrameDownLen);
    [~, tmp] = max(xFilt(idx));
    syncIdx  = tmp;  % Index relative to the subframe

    % (B) Convert correlation index to buffer index by compensating preamble length
    syncTime = round(syncIdx*ADS_B_Parameter.SyncDownsampleFactor - syncSigLen + p*subFrameLen);

    % (C) Boundary guard
    if (syncTime <= 0) || (syncTime + ADS_B_Parameter.MaxPacketLength - 1 > xLen)
      continue;
    end

    % (D) Preamble validation using chip-level energy over 16 chips
    rxSyncEnergy = energySig(syncTime + (0:syncSigLen-1));
    rxSyncSeq    = sum(reshape(rxSyncEnergy, spc, syncLen), 1);  % Energy per chip

    hi = ADS_B_Parameter.SyncSequenceHighIndices;
    lo = ADS_B_Parameter.SyncSequenceLowIndices;
    th = (mean(rxSyncSeq(hi)) + mean(rxSyncSeq(lo)))/2;          % Adaptive threshold

    % Expect chips with '1' to be > th and chips with '0' to be < th
    if all(xor((rxSyncSeq < th), ADS_B_Parameter.SyncSequence))
      % (E) Passed preamble check -> slice the 112-bit data region that follows
      packetCnt = packetCnt + 1;
      if packetCnt <= ADS_B_Parameter.MaxNumPacketsInFrame
        dataIdx = int32(ADS_B_Parameter.PreambleLength + (0:ADS_B_Parameter.LongPacketLength-1));
        packetSamples(:,packetCnt) = xBuff(syncTime + dataIdx, 1);
        syncTimeVec(packetCnt)     = syncTime;
      end
    end
  end
end

%% -------- Mode-S/ADS-B CRC-24 parity --------
function isValid = checkCRC(bits)
  % Mode S/ADS-B (DO-260) CRC-24 over 88-bit payload with 24 parity bits.
  % Generator polynomial (binary): 1111111111111010000001001  = 0x1FFF409.
  if numel(bits) < 112, isValid = false; return; end
  poly = [1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1];
  msg   = bits(1:88);
  crcRx = bits(89:112);
  dividend = [msg zeros(1,24)];
  for i = 1:numel(msg)
    if dividend(i)==1, dividend(i:i+24) = bitxor(dividend(i:i+24), poly); end
  end
  isValid = isequal(dividend(end-23:end), crcRx);
end
