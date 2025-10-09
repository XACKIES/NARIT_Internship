%% ADS-B Receiver PlutoSDR ( Model 2 : Frame Buffer + ROC + CSV + SNR )
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s)
%     -> FIR pre-filter
%     -> magnitude-squared (energy)
%     -> Frame Buffer (sliding)
%     -> matched filter (preamble template, over energy)
%     -> correlation peak search
%     -> preamble validation (chip-level energy check)
%     -> PPM demodulation (early-late)
%     -> CRC-24 parity check
%     -> ROC Curve (sample every k samples from correlation) 
%     -> CSV logging for packets (incl. SNR & AUC)            
% 
% References:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS)
%   [2] ICAO Annex 10, Vol. IV
%   [3] Mode S / ADS-B CRC-24 polynomial: 0x1FFF409
%   [4] Fawcett, "An introduction to ROC analysis", Pattern Recognition Letters, 2006
%   [5] Oppenheim & Schafer, Discrete-Time Signal Processing
%   [6] Proakis & Salehi, Digital Communications
%   [7] Huber & Ronchetti, Robust Statistics

clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT   = true;           % Enable ROC plotting
DEBUG_EVERY  = 20;             % Plot ROC every N frames
iterCount    = 0;

% CSV logging (เฉพาะ Packet CSV)
LOG_CSV   = true;
RUN_TAG   = datestr(now,'yyyymmdd_HHMMSS');
CSV_FILE  = sprintf('adsb_%s.csv', RUN_TAG);

% Create CSV headers (Packet CSV includes AUC_last)
if LOG_CSV
  [fid,msg] = fopen(CSV_FILE,'w');
  if fid == -1
    warning('adsb:csvOpen','Cannot open CSV: %s (reason: %s). Disable logging.', CSV_FILE, msg);
    LOG_CSV = false;
  else
    fprintf(fid,'timestamp,CA,ICAO_hex,SNRiq_dB,Psig_dBFS,Pnoi_dBFS,AUC\n');
    fclose(fid);
  end
end

%% ---------- ROC controls/state ----------
ROC_SAMPLE_STEP   = 100;       % sample every k samples from correlation output
ROC_POS_WINDOW_US = 2.0;       % ±2 us window around true preamble center for labels
ROC_MIN_POINTS    = 2000;      % start plotting ROC after collecting enough points
ROC_THRESH_PTS    = 200;       % number of thresholds for ROC sweep
ROC_MAX_SAMPLES   = inf;       % hard cap of collected samples
roc_scores = [];               % accumulated scores
roc_labels = [];               % accumulated labels (1/0)
AUC_last   = NaN;              % keep last computed AUC for printing/CSV

% ROC label policy when no CRC-valid this frame
ROC_LABEL_MODE = 'crc-only';   % {'crc-only','peaks-when-no-crc'}

% Pre-create ROC figure (optional)
hFigROC = [];
if DEBUG_PLOT
  hFigROC = figure(9); clf;
  set(hFigROC,'Name','ADS-B ROC','NumberTitle','off','Color','w');
  plot([0 1],[0 1],'--'); grid on; axis([0 1 0 1]); axis square;
  title('ROC — initializing...'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
  drawnow;
end

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

%% ---------- ADS-B PHY Parameters ----------
spc_chip  = round(sampRate/2e6);   % samples per 0.5 us chip
spb_bit   = 2*spc_chip;            % samples per 1.0 us bit

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

% Matched-filter template (bipolar ±1 at chip rate), time-reversed
preamble_bip = 2*SyncSequence - 1;
mf = flipud(repelem(preamble_bip(:), spc_chip));

% Sliding buffer (complex IQ)
xBuff = zeros(bufferLen,1);

conting = 0;  % count valid DF=17 across run

%% ---------- Main loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Read RF samples
  newFrame = rx(); if isempty(newFrame), continue; end

  % 2) FIR pre-filter (complex IQ)
  [newFrameFIR, firState] = filter(b_fir, 1, newFrame, firState);

  % 3) Update sliding buffer (complex IQ)
  xBuff = [xBuff(numel(newFrameFIR)+1:end); newFrameFIR];

  % 4) Instantaneous energy
  energySig = abs(xBuff).^2;

  % 5) Matched filter over energy (correlation output)
  xFilt = conv(energySig, mf, 'same');

  % 6) Packet search (returns validated packets and preamble start times)
  [packetSamples, packetCnt, syncTimeVec] = packetSearch(abs(xFilt), xBuff, energySig, ADS_B_Parameter);

  % ===== ROC accumulation each frame (either strict CRC-only or peaks-when-no-crc) =====
  pos_mask = false(size(xFilt));
  if strcmpi(ROC_LABEL_MODE,'peaks-when-no-crc')
    % ถ้าไม่มี CRC-valid ในเฟรมนี้ จะไปตั้งหน้าต่างบวกรอบพีคสูงสุดด้านล่าง
  end

  % 7) Decode each validated packet
  any_crc_valid = false;
  for k = 1:packetCnt
    region_energy = abs(packetSamples(:,k)); % energy segment (magnitude)
    nBits = floor(numel(region_energy)/spb_bit);
    if nBits < 112, continue; end

    % Early-Late PPM energy compare
    segsE = reshape(region_energy(1:nBits*spb_bit), spb_bit, nBits);
    early = sum(segsE(1:spc_chip, :), 1);
    late  = sum(segsE(spc_chip+1:end, :), 1);
    bits  = double(early > late);
    bits  = bits(1:112);

    % CRC-24 (DO-260)
    if checkCRC(bits)
      DF = bin2dec(num2str(bits(1:5)));
      if DF==17
        any_crc_valid = true;
        conting = conting+1;
        CA   = bin2dec(num2str(bits(6:8)));
        ICAO = bits(9:32);
        DATA = bits(33:88);

        % ---------- SNR (mask-only IQ around this packet) ----------
        preamble_len   = ADS_B_Parameter.PreambleLength;
        msg_len        = ADS_B_Parameter.LongPacketLength;
        syncTime       = syncTimeVec(k);  % start-of-preamble index in xBuff
        dataStart      = syncTime + preamble_len;
        dataEnd        = dataStart + msg_len - 1;
        % Boundary guard (robust)
        dataStart = max(1, dataStart);
        dataEnd   = min(length(xBuff), dataEnd);

        region_iq  = xBuff(dataStart:dataEnd); % complex IQ inside message
        P_sig_lin  = max(abs(region_iq).^2);   % peak power inside message
        % Noise: all IQ outside [syncTime - preamble_len, dataEnd]
        mask = true(size(xBuff));
        sigL = max(1, syncTime - preamble_len);
        sigR = min(length(xBuff), dataEnd);
        mask(sigL:sigR) = false;
        noise_iq = xBuff(mask);

        if isempty(noise_iq), noise_iq = xBuff; end
        P_noi_lin = mean(abs(noise_iq).^2);

        SNR_iq_dB = 10*log10(P_sig_lin/(P_noi_lin + eps));
        Psig_dBFS = 10*log10(P_sig_lin + eps);
        Pnoi_dBFS = 10*log10(P_noi_lin + eps);

        % ---------- Print pretty line (include AUC_last & Ncur) ----------
        ICAO_hex = sprintf('%06X', bin2dec(char(ICAO + '0')));
        DATA_str = char(DATA + '0');
        ts_now   = datestr(now,'HH:MM:SS.FFF');
        Ncur     = numel(roc_scores);
        fprintf(['[%s] [ADS-B] ICAO=%s | CA=%d | Message=%s | SNR=%.2f dB' ...
                 ' | Count=%d | AUC=%s | N=%d\n'], ...
                ts_now, ICAO_hex, CA, DATA_str, SNR_iq_dB, conting, ...
                ternaryStr(isnan(AUC_last),'NaN',sprintf('%.3f',AUC_last)), Ncur);

        % ---------- Packet CSV log (append) ----------
        if LOG_CSV
          tstamp = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
          [fid,msg] = fopen(CSV_FILE,'a');
          if fid ~= -1
            fprintf(fid,'%s,%d,%s,%.2f,%.2f,%.2f,%.3f\n', ...
                    tstamp, CA, ICAO_hex, SNR_iq_dB, Psig_dBFS, Pnoi_dBFS, AUC_last);
            fclose(fid);
          else
            warning('adsb:csvAppend','Packet CSV append failed: %s', msg);
          end
        end

        % ---------- Mark positive window for ROC labels (strict mode) ----------
        if strcmpi(ROC_LABEL_MODE,'crc-only')
          preamble_center = syncTime + round(preamble_len/2);
          win_samp = round((ROC_POS_WINDOW_US * 1e-6) * sampRate);
          i1 = max(1, preamble_center - win_samp);
          i2 = min(length(xFilt), preamble_center + win_samp);
          pos_mask(i1:i2) = true;
        end
      end
    end
  end

  % ---------- ROC accumulation (frame-level) ----------
  if strcmpi(ROC_LABEL_MODE,'peaks-when-no-crc') && ~any_crc_valid
    [~, loc] = max(xFilt);
    win_samp = round((ROC_POS_WINDOW_US * 1e-6) * sampRate);
    i1 = max(1, loc - win_samp); i2 = min(length(xFilt), loc + win_samp);
    pos_mask(i1:i2) = true;
  end

  % Downsample correlation to collect ROC samples
  s_idx = 1:ROC_SAMPLE_STEP:length(xFilt);
  roc_scores = [roc_scores; abs(xFilt(s_idx)).'];       %#ok<AGROW>
  roc_labels = [roc_labels; double(pos_mask(s_idx)).']; %#ok<AGROW>

  % Hard cap → finalize plotting and exit
  if numel(roc_scores) >= ROC_MAX_SAMPLES
    Ncap = ROC_MAX_SAMPLES;
    if numel(roc_scores) > Ncap
      roc_scores = roc_scores(1:Ncap);
      roc_labels = roc_labels(1:Ncap);
    end
    try
      [FPR, TPR, AUC] = computeROC(roc_scores, roc_labels, ROC_THRESH_PTS);
      AUC_last = AUC;
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
    return;
  end

  % Intermediate ROC plot (no CSV save)
  if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY)==0 && numel(roc_scores) >= ROC_MIN_POINTS
    [FPR, TPR, AUC, ~] = computeROC_verbose(roc_scores, roc_labels, ROC_THRESH_PTS);
    AUC_last = AUC;
    if isempty(hFigROC) || ~ishandle(hFigROC), hFigROC = figure(9); end
    set(0,'CurrentFigure',hFigROC); clf;
    plot(FPR, TPR, '-o', 'LineWidth', 1.25, 'MarkerSize', 3);
    grid on; axis([0 1 0 1]); axis square;
    title(sprintf('ROC (N=%d samples)  AUC = %.3f', numel(roc_scores), AUC));
    xlabel('False Positive Rate (FPR)'); ylabel('True Positive Rate (TPR)');
    drawnow limitrate;
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

%% -------- ROC compute (vectorized; returns thresholds for inspection) --------
function [FPR, TPR, AUC, thr] = computeROC_verbose(scores, labels, Nth)
  scores = scores(:);
  labels = labels(:) > 0;
  npos = sum(labels==1); nneg = sum(labels==0);
  if npos==0 || nneg==0, FPR=[0;1]; TPR=[0;1]; AUC=NaN; thr=[0;1]; return; end
  smin=min(scores); smax=max(scores);
  if smin==smax, FPR=[0;1]; TPR=[0;1]; AUC=0.5; thr=[smax;smin]; return; end

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
  thr = thr(ord);
  AUC = trapz(FPR, TPR);
end

%% -------- ROC compute (toolbox-free, minimal) --------
function [FPR, TPR, AUC] = computeROC(scores, labels, Nth)
  [FPR,TPR,AUC,~] = computeROC_verbose(scores, labels, Nth);
end

%% -------- Tiny helper --------
function s = ternaryStr(cond, a, b)
if cond, s = a; else, s = b; end
end
