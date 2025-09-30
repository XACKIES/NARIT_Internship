%% ADS-B Receiver PlutoSDR ( Model 2 : Frame Buffer + ROC )
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s) 
%     -> FIR pre-filter 
%     -> magnitude-squared (energy)
%     -> Frame Buffer (10× message window, sliding)
%     -> matched filter (preamble template)
%     -> correlation peak search 
%     -> preamble validation (chip-level energy check)
%     -> PPM demodulation (early-late)
%     -> CRC-24 parity check
%     -> ROC Curve (sample every 10 samples from correlation output)
%
% References:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS).
%   [2] ICAO Annex 10, Vol. IV.
%   [3] Mode S / ADS-B CRC-24 polynomial: 0x1FFF409.
%   [4] Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd ed.

clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT   = true;        % Enable ROC plotting
DEBUG_EVERY  = 10;          % Plot ROC every N frames
iterCount    = 0;

%% ---------- ROC controls/state ----------
ROC_SAMPLE_STEP   = 20;     % sample every 10 samples from correlation output
ROC_POS_WINDOW_US = 2.0;    % ±2 us window around true preamble center for positive labels
ROC_MIN_POINTS    = 2000;   % start plotting ROC after collecting enough points
ROC_THRESH_PTS    = 1000;    % number of thresholds for ROC sweep
ROC_MAX_SAMPLES   = 1000000; % <<< stop when total collected ROC samples reaches this cap
roc_scores = [];            % accumulated scores
roc_labels = [];            % accumulated labels (1/0)

%% ---------- Radio / Buffer ----------
fc       = 1090e6;          % ADS-B carrier frequency
sampRate = 10e6;            % 10 MS/s -> 5 samples per 0.5 us chip at 2 MHz
frameLen = 65536;           % SDR read size
bufferLen= 16*frameLen;     % Sliding buffer length

% PlutoSDR (Fast Attack AGC)
rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready (AGC Fast Attack) ... Listening 1090 MHz ...');

%% ---------- FIR front-end (pre-buffer) ----------
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286 -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];
firState = zeros(length(b_fir)-1, 1);   % streaming state

%% ---------- ADS-B PHY Parameter ----------
spc_chip  = round(sampRate/2e6);  % samples per 0.5 us chip
spb_bit   = 2*spc_chip;           % samples per 1.0 us bit

% 16-chip preamble (8 us)
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];

ADS_B_Parameter.SamplesPerChip             = spc_chip;
ADS_B_Parameter.SyncSequence               = SyncSequence;
ADS_B_Parameter.SyncSequenceLength         = length(SyncSequence);
ADS_B_Parameter.SyncSequenceHighIndices    = find(SyncSequence==1);
ADS_B_Parameter.SyncSequenceLowIndices     = find(SyncSequence==0);
ADS_B_Parameter.SyncSequenceNumHighValues  = numel(ADS_B_Parameter.SyncSequenceHighIndices);
ADS_B_Parameter.SyncSequenceNumLowValues   = numel(ADS_B_Parameter.SyncSequenceLowIndices);
ADS_B_Parameter.PreambleLength             = ADS_B_Parameter.SyncSequenceLength * spc_chip;
ADS_B_Parameter.LongPacketLength           = 112 * spb_bit;
ADS_B_Parameter.MaxPacketLength            = ADS_B_Parameter.PreambleLength + ADS_B_Parameter.LongPacketLength;
ADS_B_Parameter.MaxNumPacketsInFrame       = 64;
ADS_B_Parameter.SyncDownsampleFactor       = 1;

% Matched-filter template (bipolar ±1 at chip rate)
preamble_bip = 2*SyncSequence - 1;
mf = flipud(repelem(preamble_bip(:), spc_chip));  % time-reversed

% Sliding buffer
xBuff = zeros(bufferLen,1);

conting = 0;  % count valid DF=17

%% ---------- Main loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Read RF samples
  newFrame = rx();
  if isempty(newFrame), continue; end

  % 2) FIR pre-filter (complex IQ)
  [newFrameFIR, firState] = filter(b_fir, 1, newFrame, firState);

  % 3) Update sliding buffer
  xBuff = [xBuff(numel(newFrameFIR)+1:end); newFrameFIR];

  % 4) Instantaneous energy
  energySig = abs(xBuff).^2;

  % 5) Matched filter over energy
  xFilt = conv(energySig, mf, 'same');   % correlation output

  % 6) Packet search
  [packetSamples, packetCnt, syncTimeVec] = packetSearch(abs(xFilt), xBuff, energySig, ADS_B_Parameter);

  % 7) Decode each packet
  for k = 1:packetCnt
    region = abs(packetSamples(:,k));

    % reshape to [samples-per-bit x bits]
    nBits = floor(numel(region)/spb_bit);
    if nBits < 112, continue; end
    segs  = reshape(region(1:nBits*spb_bit), spb_bit, nBits);

    % Early-Late PPM energy compare
    early = sum(segs(1:spc_chip, :), 1);
    late  = sum(segs(spc_chip+1:end, :), 1);
    bits  = double(early > late);
    bits  = bits(1:112);

    % CRC-24 (DO-260)
    if checkCRC(bits)
      DF = bin2dec(num2str(bits(1:5)));
      if DF==17
        % ======== Your original print block (unchanged) ========
        conting = conting+1; 
        CA = bin2dec(num2str(bits(6:8))); 
        ICAO = bits(9:32); 
        DATA = bits(33:88); 
        disp('========= VALID ADS-B (DF=17) ========='); 
        fprintf('CA = %d\n', CA); 
        fprintf('ICAO = %06X\n', bin2dec(char(ICAO + '0'))); 
        fprintf('DATA = %s\n', char(DATA + '0')); 
        fprintf('#Aircraft Detecting Counting = %d\n', (conting)); 
        disp('======================================');
        % ========================================================

        % ========= ROC sampling & plotting (only ROC; no other plots) =========
        if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY)==0
          try
            % Positive window around mid-preamble of this packet
            preamble_len    = ADS_B_Parameter.PreambleLength;
            preamble_center = syncTimeVec(k) + round(preamble_len/2);

            win_samp = round((ROC_POS_WINDOW_US * 1e-6) * sampRate);
            pos_mask = false(size(xFilt));
            i1 = max(1, preamble_center - win_samp);
            i2 = min(length(xFilt), preamble_center + win_samp);
            pos_mask(i1:i2) = true;

            % Downsample correlation to collect ROC samples
            s_idx = 1:ROC_SAMPLE_STEP:length(xFilt);
            roc_scores = [roc_scores; abs(xFilt(s_idx)).'];      %#ok<AGROW>
            roc_labels = [roc_labels; double(pos_mask(s_idx)).']; %#ok<AGROW>

            % ---- HARD CAP: finalize when reaching ROC_MAX_SAMPLES ----
            if numel(roc_scores) >= ROC_MAX_SAMPLES
              Ncap = ROC_MAX_SAMPLES;
              if numel(roc_scores) > Ncap
                roc_scores = roc_scores(1:Ncap);
                roc_labels = roc_labels(1:Ncap);
              end
              % Final ROC plot with exact N = ROC_MAX_SAMPLES
              try
                [FPR, TPR, AUC] = computeROC(roc_scores, roc_labels, ROC_THRESH_PTS);
                figure(9); clf;
                plot(FPR, TPR, '-o', 'LineWidth', 1.25, 'MarkerSize', 3);
                grid on; axis([0 1 0 1]); axis square;
                title(sprintf('ROC (N=%d samples)  AUC = %.3f', Ncap, AUC));
                xlabel('False Positive Rate (FPR)'); ylabel('True Positive Rate (TPR)');
                drawnow;
              catch
              end
              fprintf('🛑 Reached ROC_MAX_SAMPLES = %d. Stopping...\n', Ncap);
              try, release(rx); catch, end
              return; % exit after final plot
            end
            % ---- Intermediate ROC plot (optional) ----
            if numel(roc_scores) >= ROC_MIN_POINTS
              [FPR, TPR, AUC] = computeROC(roc_scores, roc_labels, ROC_THRESH_PTS);
              figure(9); clf;
              plot(FPR, TPR, '-o', 'LineWidth', 1.25, 'MarkerSize', 3);
              grid on; axis([0 1 0 1]); axis square;
              title(sprintf('ROC (N=%d samples)  AUC = %.3f', numel(roc_scores), AUC));
              xlabel('False Positive Rate (FPR)'); ylabel('True Positive Rate (TPR)');
              drawnow limitrate;
            end
          catch
          end
        end
        % =====================================================================
      end
    end
  end
end

%% -------- helper: Framework-like packet search --------
function [packetSamples, packetCnt, syncTimeVec] = packetSearch(xFilt, xBuff, energySig, ADS_B_Parameter)
  spc = ADS_B_Parameter.SamplesPerChip;
  syncLen    = ADS_B_Parameter.SyncSequenceLength;
  syncSigLen = syncLen*spc;
  xLen       = length(xBuff);

  subFrameLen     = ADS_B_Parameter.MaxPacketLength;
  subFrameDownLen = subFrameLen / ADS_B_Parameter.SyncDownsampleFactor;
  numSubFrames    = floor(xLen / subFrameLen);

  packetSamples = zeros(ADS_B_Parameter.LongPacketLength, ADS_B_Parameter.MaxNumPacketsInFrame, 'like', xBuff);
  syncTimeVec   = zeros(ADS_B_Parameter.MaxNumPacketsInFrame,1);
  packetCnt     = 0;

  for p = 0:(numSubFrames-2)
    % (A) Peak search inside this subframe
    idx = double(p)*subFrameDownLen + (1:subFrameDownLen);
    [~, tmp] = max(xFilt(idx));
    syncIdx  = tmp;  % index within subframe

    % (B) Convert to buffer index; compensate preamble length
    syncTime = round(syncIdx*ADS_B_Parameter.SyncDownsampleFactor - syncSigLen + p*subFrameLen);

    % (C) Boundary guard
    if (syncTime <= 0) || (syncTime + ADS_B_Parameter.MaxPacketLength - 1 > xLen)
      continue;
    end

    % (D) Preamble validation at chip-level (energy)
    rxSyncEnergy = energySig(syncTime + (0:syncSigLen-1));
    rxSyncSeq    = sum(reshape(rxSyncEnergy, spc, syncLen), 1);

    hi = ADS_B_Parameter.SyncSequenceHighIndices;
    lo = ADS_B_Parameter.SyncSequenceLowIndices;
    th = (mean(rxSyncSeq(hi)) + mean(rxSyncSeq(lo)))/2;

    if all(xor((rxSyncSeq < th), ADS_B_Parameter.SyncSequence))
      % (E) Valid → slice following 112-bit data
      packetCnt = packetCnt + 1;
      if packetCnt <= ADS_B_Parameter.MaxNumPacketsInFrame
        dataIdx = int32(ADS_B_Parameter.PreambleLength + (0:ADS_B_Parameter.LongPacketLength-1));
        packetSamples(:,packetCnt) = xBuff(syncTime + dataIdx, 1);
        syncTimeVec(packetCnt)     = syncTime;  % start-of-preamble in buffer
      end
    end
  end
end

%% -------- Mode-S/ADS-B CRC-24 parity --------
function isValid = checkCRC(bits)
  % CRC over 88-bit payload with 24 parity (DO-260), poly = 0x1FFF409
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

%% -------- ROC compute (toolbox-free) --------
function [FPR, TPR, AUC] = computeROC(scores, labels, Nth)
  scores = scores(:);
  labels = labels(:) > 0;
  npos = sum(labels==1); nneg = sum(labels==0);
  if npos==0 || nneg==0, FPR=[0;1]; TPR=[0;1]; AUC=NaN; return; end
  smin=min(scores); smax=max(scores);
  if smin==smax, FPR=[0;1]; TPR=[0;1]; AUC=0.5; return; end

  thr = linspace(smax, smin, max(3, Nth)).';
  TPR = zeros(numel(thr),1); FPR = zeros(numel(thr),1);
  for i=1:numel(thr)
    t  = thr(i);
    yp = (scores >= t);
    TP = sum( yp &  labels );
    FP = sum( yp & ~labels );
    FN = sum(~yp &  labels );
    TN = sum(~yp & ~labels );
    TPR(i) = TP / (TP + FN + eps);
    FPR(i) = FP / (FP + TN + eps);
  end
  [FPR, ord] = sort(FPR);
  TPR = TPR(ord);
  AUC = trapz(FPR, TPR);
end
