% ADS-B Receiver PlutoSDR (Model 1: Direct Xcorr + FIR + Adaptive Threshold)
% + Proper Alignment + ROC (plot-only) + CSV Logging + No logical-scalar errors
%
% Chain:
%   PlutoSDR RX (10 MS/s)
%   -> complex FIR front-end (linear-phase)
%   -> magnitude
%   -> matched filter (preamble template, direct xcorr)
%   -> adaptive threshold (median + k*MAD)
%   -> peak pick
%   -> ALIGN (+L/2 from correlation peak)
%   -> PPM demod
%   -> CRC-24 parity check
%   -> Packet CSV log (includes AUC_last)
%   -> ROC (scores sampled, plot-only; print AUC)
% 
% Refs: [1] DO-260B/C, [2] ICAO Annex10 Vol.IV, [3] CRC24 0x1FFF409,
%       [4] Oppenheim & Schafer, [5] Proakis & Salehi, [6] Huber & Ronchetti

clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT      = true;          % show ROC figure
DEBUG_EVERY     = 10;            % plot every N frames
th_k            = 5;             % threshold = median + k*MAD

% CSV logging (Packet only)
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
ROC_SAMPLE_STEP   = 100;          % sampling step on correlation output
ROC_POS_WINDOW_US = 2.0;          % ± window around centers
ROC_MIN_POINTS    = 200;          % show curve sooner
ROC_THRESH_PTS    = 200;          % thresholds per sweep
ROC_MAX_SAMPLES   = inf;          % cap
roc_scores = [];
roc_labels = [];
AUC_last   = NaN;                 % keep last computed AUC for printing/CSV

% Label policy when no CRC-valid this frame:
ROC_LABEL_MODE = 'peaks-when-no-crc';  % {'crc-only','peaks-when-no-crc'}

iterCount = 0;

% Pre-create ROC figure handle (avoid figure focus issues)
hFigROC = [];
if DEBUG_PLOT
  hFigROC = figure(9); clf;
  set(hFigROC,'Name','ADS-B ROC','NumberTitle','off','Color','w');
  plot([0 1],[0 1],'--'); grid on; axis([0 1 0 1]); axis square;
  title('ROC — initializing...'); xlabel('False Positive Rate'); ylabel('True Positive Rate');
  drawnow;
end

%% ---------- Radio / Buffer ----------
fc        = 1090e6;
sampRate  = 10e6;               % 10 MS/s
frameLen  = 65536;              % ~6.55 ms

rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');
conting = 0;

%% ---------- FIR front-end (linear-phase) ----------
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286  -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];
zi_fir = zeros(numel(b_fir)-1,1);

%% ---------- ADS-B PHY Parameters ----------
os = round(sampRate/1e6);               % ~10 samp/us
SyncSequence   = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];
SamplesPerChip = max(1, round(os/2));   % 0.5 us per chip
preamble       = repelem(SyncSequence, SamplesPerChip);
mf             = flipud(preamble(:));
L              = numel(preamble);

half = floor(os/2);
bit0 = [ones(1,half) zeros(1,os-half)];
bit1 = [zeros(1,half) ones(1,os-half)];
msgBits = 112;

%% ---------- Main Loop ----------
while true
  iterCount = iterCount + 1;
  rxSig = rx(); if isempty(rxSig); continue; end

  % 1) FIR (complex)
  [y_fir, zi_fir] = filter(b_fir, 1, rxSig, zi_fir);

  % 2) Matched filter on magnitude
  c = abs(conv(abs(y_fir), mf, 'same'));

  % 3) Robust threshold
  med = median(c);
  madv = median(abs(c - med)) + eps;
  th  = med + th_k * madv;

  % 4) Detections
  locs = find(c > th);
  true_peaks_this_frame = [];

  if ~isempty(locs)
    for st = locs(:).'                     % ensure scalar 'st'
      % 5) ALIGN start
      st_start   = st + floor(L/2);
      msgSamples = msgBits * os;

      % --- Harden bound check (vector-safe) ---
      if any( (st_start < 1) | ((st_start + msgSamples - 1) > numel(y_fir)) )
        continue;
      end

      % 6) PPM demod (magnitude)
      region_mag = abs(y_fir(st_start : st_start + msgSamples - 1));
      segs   = reshape(region_mag, os, msgBits);
      score0 = sum(segs .* bit0.', 1);
      score1 = sum(segs .* bit1.', 1);
      bits   = double(score1 > score0);

      % 7) CRC check
      DF = bin2dec(char(bits(1:5) + '0'));
      if (DF == 17) && checkCRC(bits)
        CA   = bin2dec(char(bits(6:8)    + '0'));
        ICAO = bits(9:32);
        DATA = bits(33:88);
        conting = conting + 1;

        % --- SNR (mask-only, IQ) ---
        region_iq  = y_fir(st_start : st_start + msgSamples - 1);
        P_sig_lin  = max(abs(region_iq).^2);
        mask = true(size(y_fir));
        sigL = max(1, st - L);
        sigR = min(numel(y_fir), st_start + msgSamples);
        mask(sigL:sigR) = false;
        noise_iq = y_fir(mask);
        if isempty(noise_iq), noise_iq = y_fir; end
        P_noi_lin = mean(abs(noise_iq).^2);
        SNR_iq_dB = 10*log10(P_sig_lin/(P_noi_lin + eps));
        Psig_dBFS = 10*log10(P_sig_lin + eps);
        Pnoi_dBFS = 10*log10(P_noi_lin + eps);

        % --- Pretty packet log line ---
        Ncur     = numel(roc_scores);  % << define before use
        ICAO_hex = sprintf('%06X', bin2dec(char(ICAO + '0')));
        DATA_str = char(DATA + '0');
        ts_now   = datestr(now,'HH:MM:SS.FFF');
        fprintf('[%s] [ADS-B] ICAO=%s | CA=%d | Message=%s | SNR=%.2f dB | Count=%d | AUC=%s  (N=%d)\n' , ...
                ts_now, ICAO_hex, CA, DATA_str, SNR_iq_dB, conting, ...
                ternaryStr(isnan(AUC_last),'NaN',sprintf('%.3f',AUC_last)) , Ncur);

        % --- Packet CSV log (append AUC_last) ---
        if LOG_CSV
          ts = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
          [fid,msg] = fopen(CSV_FILE,'a');
          if fid ~= -1
            fprintf(fid,'%s,%d,%s,%.2f,%.2f,%.2f,%.3f\n', ...
                    ts, CA, ICAO_hex, SNR_iq_dB, Psig_dBFS, Pnoi_dBFS, AUC_last);
            fclose(fid);
          else
            warning('adsb:csvAppend','Packet CSV append failed: %s', msg);
          end
        end

        % mark validated center
        true_peaks_this_frame(end+1) = st; %#ok<AGROW>
      end
    end
  end

  %% ---------- ROC accumulation (every frame; vector-safe) ----------
  pos_mask = false(size(c));
  if ~isempty(true_peaks_this_frame)
    win_samp = round((ROC_POS_WINDOW_US * 1e-6) * sampRate);
    for k = 1:numel(true_peaks_this_frame)
      center = true_peaks_this_frame(k);
      i1 = max(1, center - win_samp);
      i2 = min(numel(c), center + win_samp);
      pos_mask(i1:i2) = true;
    end
  elseif strcmpi(ROC_LABEL_MODE,'peaks-when-no-crc') && ~isempty(locs)
    win_samp = round((ROC_POS_WINDOW_US * 1e-6) * sampRate);
    for k = 1:numel(locs)
      center = locs(k);
      i1 = max(1, center - win_samp);
      i2 = min(numel(c), center + win_samp);
      pos_mask(i1:i2) = true;
    end
  end

  s_idx = 1:ROC_SAMPLE_STEP:numel(c);
  roc_scores = [roc_scores; c(s_idx).'];
  roc_labels = [roc_labels; double(pos_mask(s_idx)).'];

  %% ---------- ROC plot-only (no CSV save) ----------
  if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY) == 0
    Ncur = numel(roc_scores);
    npos = sum(roc_labels==1);
    nneg = sum(roc_labels==0);

    cond_full = (Ncur >= ROC_MIN_POINTS) && (npos>0) && (nneg>0);
    Nth_use = cond_full ? ROC_THRESH_PTS : max(10, round(ROC_THRESH_PTS/5));

    [FPR, TPR, AUC, ~] = computeROC(roc_scores, roc_labels, Nth_use);
    AUC_last = AUC;

    if isempty(hFigROC) || ~ishandle(hFigROC), hFigROC = figure(9); end
    set(0,'CurrentFigure',hFigROC); clf;
    plot(FPR, TPR, '-o', 'LineWidth', 1.25, 'MarkerSize', 3);
    grid on; axis([0 1 0 1]); axis square;
    if isnan(AUC)
      title(sprintf('ROC (N=%d) — collecting labels...', Ncur));
    else
      title(sprintf('ROC (N=%d)  AUC = %.3f', Ncur, AUC));
    end
    xlabel('False Positive Rate'); ylabel('True Positive Rate');
    drawnow;

    if Ncur >= ROC_MAX_SAMPLES
      try, release(rx); catch, end
      return;
    end
  end
end

%% ---------- Helper: Mode-S/ADS-B CRC-24 ----------
function isValid = checkCRC(bits)
  % Polynomial: 0x1FFF409
  poly     = [1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1];
  msg      = bits(1:88);
  crcRx    = bits(89:112);
  dividend = [msg zeros(1,24)];
  for i = 1:length(msg)
    if dividend(i) == 1
      dividend(i:i+24) = bitxor(dividend(i:i+24), poly);
    end
  end
  crcCalc = dividend(end-23:end);
  isValid = isequal(crcCalc, crcRx);
end

%% ---------- Helper: ROC computation ----------
function [FPR, TPR, AUC, thresh] = computeROC(scores, labels, Nth)
  scores = scores(:);
  labels = labels(:) > 0;
  npos = sum(labels==1); nneg = sum(labels==0);
  if (npos==0) || (nneg==0), FPR=[0;1]; TPR=[0;1]; AUC=NaN; thresh=[0;0]; return; end
  smin=min(scores); smax=max(scores);
  if smin==smax, FPR=[0;1]; TPR=[0;1]; AUC=0.5; thresh=[0;0]; return; end
  Nth = max(3, Nth);
  thresh=linspace(smax,smin,Nth).';
  TPR=zeros(numel(thresh),1); FPR=zeros(numel(thresh),1);
  for i=1:numel(thresh)
    t=thresh(i);
    yp=(scores>=t);
    TP=sum( yp &  labels );
    FP=sum( yp & ~labels );
    FN=sum(~yp &  labels );
    TN=sum(~yp & ~labels );
    TPR(i)=TP/(TP+FN+eps);
    FPR(i)=FP/(FP+TN+eps);
  end
  [FPR,order]=sort(FPR);
  TPR=TPR(order);
  thresh=thresh(order);
  AUC=trapz(FPR,TPR);
end

%% ---------- Tiny helper ----------
function s = ternaryStr(cond, a, b)
% ternary string helper: if cond is true -> a, else -> b
if cond, s = a; else, s = b; end
end
