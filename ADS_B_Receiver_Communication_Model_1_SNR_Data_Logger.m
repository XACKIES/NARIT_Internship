%% ADS-B Receiver PlutoSDR (Model 1: Direct Xcorr + FIR + Adaptive Threshold + Proper Alignment)
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s)
%     -> complex FIR front-end (linear-phase, streaming state via z_i/z_f)
%     -> magnitude (abs)
%     -> matched filter (preamble template, direct xcorr via conv)
%     -> adaptive threshold (median + k*MAD)
%     -> peak pick (candidates above threshold)
%     -> slice after preamble (ALIGN: +L/2 from correlation peak)
%     -> PPM demodulation (bit0/bit1 templates)
%     -> CRC-24 parity check
%
% Engineering references:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS).
%   [2] ICAO Annex 10, Vol. IV.
%   [3] Mode S / ADS-B CRC-24 polynomial: 0x1FFF409 (DO-260).
%   [4] Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd ed.
%   [5] Proakis & Salehi, Digital Communications, 5th ed.
%   [6] Kay, Estimation Theory, Prentice Hall.
%
% Notes:
%   - Requires: Communications Toolbox Support Package for ADALM-PLUTO Radio.
%   - AGC Fast Attack recommended at 10 MS/s for bursty ADS-B.

clear; clc;

%% ---------- User Controls ----------
iterCount    = 0;
th_k         = 5.0;       % Adaptive threshold parameter (median + k*MAD)

%% ---------- CSV Logging ----------
LOG_CSV   = true;
RUN_TAG   = datestr(now,'yyyymmdd_HHMMSS');
CSV_FILE  = sprintf('adsb_%s.csv', RUN_TAG);
if LOG_CSV
  [fid,msg] = fopen(CSV_FILE,'w');
  if fid == -1
    warning('Cannot open CSV for write: %s (reason: %s). Logging disabled.', CSV_FILE, msg);
    LOG_CSV = false;
  else
    fprintf(fid,'timestamp,CA,ICAO_hex,SNRiq_dB,Psig_dBFS,Pnoi_dBFS\n');
    fclose(fid);
  end
end

%% ---------- Radio / Buffer ----------
fc        = 1090e6;        % ADS-B carrier
sampRate  = 10e6;          % 10 MS/s
frameLen  = 65536;         % Samples per frame

rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');
conting = 0;

%% ---------- FIR front-end ----------
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286 -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];
M_fir  = numel(b_fir);
zi_fir = zeros(M_fir-1,1);

%% ---------- ADS-B PHY Parameters ----------
os = round(sampRate/1e6);
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];
SamplesPerChip = max(1, round(os/2));
preamble       = repelem(SyncSequence, SamplesPerChip);
mf             = flipud(preamble(:));
half           = floor(os/2);
bit0           = [ones(1,half) zeros(1,os-half)];
bit1           = [zeros(1,half) ones(1,os-half)];
L              = numel(preamble);
msgBits        = 112;

%% ---------- Main Loop ----------
while true
  iterCount = iterCount + 1;

  rxSig = rx();
  if isempty(rxSig), continue; end

  [y_fir, zi_fir] = filter(b_fir, 1, rxSig, zi_fir);
  c = abs(conv(abs(y_fir), mf, 'same'));

  med = median(c);
  mad = median(abs(c - med)) + eps;
  th  = med + th_k * mad;

  locs = find(c > th);
  if isempty(locs), continue; end

  for st = locs(:)'
    st_start   = st + floor(L/2);
    msgSamples = msgBits * os;
    if st_start < 1 || st_start + msgSamples - 1 > length(y_fir), continue; end

    region_mag = abs(y_fir(st_start : st_start + msgSamples - 1));
    region_iq  =       y_fir(st_start : st_start + msgSamples - 1);

    segs   = reshape(region_mag, os, msgBits);
    score0 = sum(segs .* bit0.', 1);
    score1 = sum(segs .* bit1.', 1);
    bits   = double(score1 > score0);

    DF = bin2dec(char(bits(1:5) + '0'));
    if DF == 17 && checkCRC(bits)
      CA   = bin2dec(char(bits(6:8)    + '0'));
      ICAO = bits(9:32);
      DATA = bits(33:88);
      conting = conting+1;

      % ---------- SNR (mask-only) ----------
      P_sig_lin = max(abs(region_iq).^2);
      mask = true(size(y_fir));
      sigL = max(1, st - L);
      sigR = min(length(y_fir), st_start + msgSamples);
      mask(sigL:sigR) = false;
      noise_iq = y_fir(mask);
      if isempty(noise_iq), noise_iq = y_fir; end
      P_noi_lin = mean(abs(noise_iq).^2);
      SNR_iq_dB = 10*log10(P_sig_lin/(P_noi_lin + eps));
      Psig_dBFS = 10*log10(P_sig_lin + eps);
      Pnoi_dBFS = 10*log10(P_noi_lin + eps);

      % ---------- Print ----------
      disp('========= VALID ADS-B (DF=17) =========');
      fprintf('CA   = %d\n', CA);
      fprintf('ICAO = %06X\n', bin2dec(char(ICAO + '0')));
      fprintf('DATA = %s\n', char(DATA + '0'));
      fprintf('#Aircraft Detecting Counting = %d\n', conting);
      fprintf('SNR (IQ, mask-only) = %.2f dB\n', SNR_iq_dB);
      %fprintf('Psig ~ %.2f dBFS, Pnoi ~ %.2f dBFS\n', Psig_dBFS, Pnoi_dBFS);
      disp('======================================');

      % ---------- Log to CSV ----------
      if LOG_CSV
        ts = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
        ICAO_hex = sprintf('%06X', bin2dec(char(ICAO + '0')));
        [fid,msg] = fopen(CSV_FILE,'a');
        if fid ~= -1
          fprintf(fid,'%s,%d,%s,%.2f,%.2f,%.2f\n', ts, CA, ICAO_hex, SNR_iq_dB, Psig_dBFS, Pnoi_dBFS);
          fclose(fid);
        else
          warning('CSV append failed: %s', msg);
        end
      end
    end
  end
end

%% ---------- Helper: Mode-S/ADS-B CRC-24 ----------
function isValid = checkCRC(bits)
  poly = [1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1];
  msg  = bits(1:88);
  crcRx = bits(89:112);
  dividend = [msg zeros(1,24)];
  for i = 1:length(msg)
    if dividend(i) == 1
      dividend(i:i+24) = bitxor(dividend(i:i+24), poly);
    end
  end
  crcCalc = dividend(end-23:end);
  isValid = isequal(crcCalc, crcRx);
end
