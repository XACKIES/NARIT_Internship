%% ADS-B Receiver PlutoSDR (Model 1: Direct Xcorr + FIR + Adaptive Threshold + ROC only under DEBUG)
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s)
%     -> FIR front-end
%     -> magnitude
%     -> matched filter (preamble template, direct xcorr)
%     -> adaptive threshold (median + k*MAD)
%     -> peak pick
%     -> PPM demod
%     -> CRC-24 parity check
%     -> ROC Curve (sampling every N samples, plotted only under DEBUG)
%
% References:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS).
%   [2] ICAO Annex 10, Vol. IV.
%   [3] DO-260 Mode S / ADS-B CRC-24 polynomial: 0x1FFF409.
%   [4] Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd ed.
%   [5] Huber, Robust Statistics (Median/MAD thresholding).

clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT      = true;       % Enable/disable ROC plotting
DEBUG_EVERY     = 10;         % Plot every N frames
iterCount       = 0;
th_k            = 6.0;        % Median + k*MAD thresholding factor

%% ---------- ROC controls/state ----------
ROC_SAMPLE_STEP   = 100;       % Collect correlation scores every [ROC_SAMPLE_STEP] samples
ROC_POS_WINDOW_US = 2.0;      % Positive label window ±2 us around preamble center
ROC_MIN_POINTS    = 2000;     % Minimum points before plotting ROC
ROC_THRESH_PTS    = 1000;      % Number of thresholds for ROC curve sweep
ROC_MAX_SAMPLES   = 1000000;   % Maximum samples to collect for ROC
roc_scores = [];              % Collected correlation magnitudes
roc_labels = [];              % Collected labels (1=positive, 0=negative)

%% ---------- Radio / Buffer ----------
fc        = 1090e6;           % ADS-B center frequency
sampRate  = 10e6;             % 10 MS/s
frameLen  = 65536;            % Frame size

rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');
conting = 0;                  % Aircraft detection counter

%% ---------- FIR front-end ----------
% Linear-phase FIR (symmetric real coefficients)
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286 -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];
M_fir  = numel(b_fir);
zi_fir = zeros(M_fir-1,1);     % Streaming filter state

%% ---------- ADS-B PHY Parameters ----------
os = round(sampRate/1e6);      % Oversampling factor (~10 samples per µs)
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];
SamplesPerChip = max(1, round(os/2)); % 0.5 µs chip
preamble = repelem(SyncSequence, SamplesPerChip);
mf       = flipud(preamble(:));       % Matched filter kernel
L        = numel(preamble);

half = floor(os/2);
bit0 = [ones(1,half) zeros(1,os-half)]; % Template for bit-0
bit1 = [zeros(1,half) ones(1,os-half)]; % Template for bit-1

%% ---------- Main Loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Receive RF samples
  rxSig = rx();
  if isempty(rxSig), continue; end

  % 2) Apply FIR filter (complex in/out)
  [y_fir, zi_fir] = filter(b_fir, 1, rxSig, zi_fir);

  % 3) Matched filter on magnitude
  c = abs(conv(abs(y_fir), mf, 'same'));

  % 4) Adaptive threshold (median + k*MAD)
  med = median(c);
  mad = median(abs(c - med)) + eps;
  th  = med + th_k * mad;

  % 5) Candidate detections above threshold
  locs = find(c > th);

  true_peaks_this_frame = [];
  if ~isempty(locs)
    for st = locs(:)'
      % 6) Align payload start = st + L/2
      st_start = st + floor(L/2);

      % 7) Guard against buffer overflow
      msgSamples = 112 * os;
      if st_start < 1 || st_start + msgSamples - 1 > length(y_fir), continue; end

      % 8) PPM demodulation
      region = abs(y_fir(st_start : st_start + msgSamples - 1));
      segs   = reshape(region, os, 112);
      score0 = sum(segs .* bit0.', 1);
      score1 = sum(segs .* bit1.', 1);
      bits   = double(score1 > score0);

      % 9) CRC check
      DF = bin2dec(char(bits(1:5) + '0'));
      if DF == 17 && checkCRC(bits)
        CA   = bin2dec(char(bits(6:8) + '0'));
        ICAO = bits(9:32);
        DATA = bits(33:88);
        conting = conting + 1;

        disp('========= VALID ADS-B (DF=17) =========');
        fprintf('CA   = %d\n', CA);
        fprintf('ICAO = %06X\n', bin2dec(char(ICAO + '0')));
        fprintf('DATA = %s\n', char(DATA + '0'));
        fprintf('#Aircraft Detecting Counting = %d\n', conting);
        disp('======================================');

        true_peaks_this_frame(end+1) = st; %#ok<AGROW>

        % ---------- ROC accumulation (only under DEBUG) ----------
        if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY) == 0
          try
            % (A) Positive label mask around true preamble centers
            win_samp = round((ROC_POS_WINDOW_US * 1e-6) * sampRate);
            pos_mask = false(size(c));
            for kkk = 1:numel(true_peaks_this_frame)
              center = true_peaks_this_frame(kkk);
              i1 = max(1, center - win_samp);
              i2 = min(length(c), center + win_samp);
              pos_mask(i1:i2) = true;
            end

            % (B) Downsample correlation values for ROC
            s_idx = 1:ROC_SAMPLE_STEP:length(c);
            roc_scores = [roc_scores; c(s_idx).'];                  %#ok<AGROW>
            roc_labels = [roc_labels; double(pos_mask(s_idx)).'];   %#ok<AGROW>

            % (C) If maximum reached: trim, plot final ROC, and stop
            if numel(roc_scores) >= ROC_MAX_SAMPLES
              Ncap = ROC_MAX_SAMPLES;
              if numel(roc_scores) > Ncap
                roc_scores = roc_scores(1:Ncap);
                roc_labels = roc_labels(1:Ncap);
              end
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
              return;
            end

            % (D) Intermediate plotting once enough points are collected
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
      end
    end
  end
end

%% ---------- Helper: CRC-24 ----------
function isValid = checkCRC(bits)
% Mode S / ADS-B (DO-260) CRC-24 over 88-bit payload + 24 parity bits
% Generator polynomial: 0x1FFF409
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
function [FPR, TPR, AUC] = computeROC(scores, labels, Nth)
% Compute ROC curve given scores and binary labels
scores = scores(:);
labels = labels(:) > 0;
npos = sum(labels==1); nneg = sum(labels==0);
if npos==0 || nneg==0
  FPR=[0;1]; TPR=[0;1]; AUC=NaN; return;
end
smin=min(scores); smax=max(scores);
if smin==smax
  FPR=[0;1]; TPR=[0;1]; AUC=0.5; return;
end
thresh=linspace(smax,smin,max(3,Nth)).';
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
AUC=trapz(FPR,TPR);
end
