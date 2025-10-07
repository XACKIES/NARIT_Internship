%% ADS-B Receiver PlutoSDR (Direct Xcorr + FIR + MAD Threshold + Proper Alignment)
% Local CPR (robust, reference-only), SNR (mask-only IQ) [ONLY WHEN POS+RANGE],
% CSV logging, great-circle range from a fixed reference.
% NOTE: The main PHY chain is unchanged. Positioning/SNR/CSV are appended.

clear; clc;

%% ---------- User Controls ----------
iterCount = 0;
th_k      = 5.0;   % detection threshold = median + k*MAD

% Reference site (Chiang Mai). Set to your actual antenna coordinates (WGS-84).
REF_LAT = 18.851836;
REF_LON = 98.957535;

% Discard position solutions farther than this range from the reference (km).
MAX_RANGE_KM = 500;

%% ---------- CSV Logging (timestamped filename) ----------
LOG_CSV  = true;
RUN_TAG  = datestr(now,'yyyymmdd_HHMMSS');                   % e.g., 20251007_154215
CSV_FILE = sprintf('adsb_%s.csv', RUN_TAG);
if LOG_CSV
  [fid,msg] = fopen(CSV_FILE,'w');
  if fid == -1
    warning('Cannot open CSV for write: %s (reason: %s). Logging disabled.', CSV_FILE, msg);
    LOG_CSV = false;
  else
    % Only rows with valid position+range are written
    fprintf(fid,'timestamp,ICAO_hex,DF,TC,SNRiq_dB,Psig_dBFS,Pnoi_dBFS,lat_deg,lon_deg,range_km,src\n');
    fclose(fid);
  end
end

%% ---------- Radio / Buffer ----------
fc        = 1090e6;   % ADS-B carrier
sampRate  = 10e6;     % 10 MS/s
frameLen  = 65536;    % samples per frame

rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');
conting = 0;

%% ---------- FIR front-end ----------
% Linear-phase FIR (symmetric). Streaming with internal state.
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286 -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];
M_fir  = numel(b_fir);
zi_fir = zeros(M_fir-1,1);

%% ---------- ADS-B PHY Parameters ----------
os = round(sampRate/1e6);            % ~10 samples per 1 µs bit
SyncSequence   = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];
SamplesPerChip = max(1, round(os/2));     % 0.5 µs per chip
preamble       = repelem(SyncSequence, SamplesPerChip);
mf             = flipud(preamble(:));     % matched-filter kernel
half           = floor(os/2);
bit0           = [ones(1,half) zeros(1,os-half)]; % PPM: pulse in first half
bit1           = [zeros(1,half) ones(1,os-half)]; % PPM: pulse in second half
L              = numel(preamble);          % preamble length in samples
msgBits        = 112;                      % DF17 bits (incl. parity)

%% ---------- Main Loop ----------
while true
  

  % 1) Acquire one frame from Pluto
  rxSig = rx();                     if isempty(rxSig), continue; end

  % 2) FIR (streaming state)
  [y_fir, zi_fir] = filter(b_fir, 1, rxSig, zi_fir);

  % 3) Matched filter on magnitude to find Mode S preamble
  c = abs(conv(abs(y_fir), mf, 'same'));

  % 4) Robust threshold: median + k·MAD
  med = median(c);
  mad = median(abs(c - med)) + eps;
  th  = med + th_k * mad;

  % 5) Candidate peaks above threshold
  locs = find(c > th);              if isempty(locs), continue; end

  % 6) For each candidate
  for st = locs(:)'   % st ≈ center of preamble (conv 'same' → group delay ≈ L/2)
    st_start   = st + floor(L/2);                 % start of payload (after preamble)
    msgSamples = msgBits * os;
    if st_start < 1 || st_start + msgSamples - 1 > length(y_fir), continue; end

    % 7) Extract payload window (magnitude for demod; IQ for SNR)
    region_mag = abs(y_fir(st_start : st_start + msgSamples - 1));
    region_iq  =       y_fir(st_start : st_start + msgSamples - 1);

    % 8) PPM demodulation via template scores
    segs   = reshape(region_mag, os, msgBits);
    score0 = sum(segs .* bit0.', 1);
    score1 = sum(segs .* bit1.', 1);
    bits   = double(score1 > score0);    % bit=1 if second half is stronger

    % 9) CRC check → accept DF=17 only
    DF = bin2dec(char(bits(1:5) + '0'));
    if ~(DF == 17 && checkCRC(bits)), continue; end

    conting = conting + 1;
    ICAO     = bits(9:32);
    ICAO_hex = sprintf('%06X', bin2dec(char(ICAO + '0')));
    ME       = bits(33:88);
    TC       = bin2dec(char(ME(1:5) + '0'));

    % === Local CPR (robust, single-frame) for Airborne Position (TC 9–18) ===
    lat_deg = NaN; lon_deg = NaN; range_km = NaN; src = "local-robust";
    if TC >= 9 && TC <= 18
      % ME fields: [1..5]=TC, [6..8]=SS, [9..20]=Alt/Q, [21]=T, [22]=F,
      %            [23..39]=LatCPR, [40..56]=LonCPR (17 bits each)
      F      = ME(22);                                  % 0 even, 1 odd frame
      latCPR = bin2dec(char(ME(23:39) + '0'));
      lonCPR = bin2dec(char(ME(40:56) + '0'));
      Y = latCPR/131072;                                % 2^17
      X = lonCPR/131072;

      [okL, lg, ln] = cprDecodeLocalRobustSimple(REF_LAT, REF_LON, Y, X, F, MAX_RANGE_KM);
      if okL
        lat_deg  = lg;
        lon_deg  = ln;
        range_km = haversine_km(REF_LAT, REF_LON, lat_deg, lon_deg);
      end
    end

    % === SNR (mask-only IQ) — computed always, but reported/logged only with POS+RANGE ===
    % Signal power: peak(|IQ|^2) within the payload window (robust to timing jitter)
    P_sig_lin = max(abs(region_iq).^2);
    % Noise power: mean(|IQ|^2) outside [preamble..payload] in the same frame
    mask = true(size(y_fir));
    sigL = max(1,               st - L);
    sigR = min(length(y_fir),   st_start + msgSamples);
    mask(sigL:sigR) = false;
    noise_iq = y_fir(mask);     if isempty(noise_iq), noise_iq = y_fir; end
    P_noi_lin = mean(abs(noise_iq).^2);

    SNR_iq_dB = 10*log10(P_sig_lin/(P_noi_lin + eps));
    Psig_dBFS = 10*log10(P_sig_lin + eps);
    Pnoi_dBFS = 10*log10(P_noi_lin + eps);

    % === Report & CSV only when position+range are valid ===
    ts = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS'));
    if ~isnan(lat_deg) && ~isnan(range_km)
      iterCount = iterCount + 1;  
      disp('========= VALID ADS-B (DF=17) =========');
      fprintf('Time = %s\n', ts);
      fprintf('ICAO = %s  DF=17  TC=%d\n', ICAO_hex, TC);
      fprintf('POS(lat,lon) = %.6f°, %.6f°  (Local CPR @ ref %.6f, %.6f)  Range ≈ %.3f km\n', ...
              lat_deg, lon_deg, REF_LAT, REF_LON, range_km);
      fprintf('SNR (IQ, mask-only) = %.2f dB  | Psig=%.2f dBFS  Pnoi=%.2f dBFS\n', SNR_iq_dB, Psig_dBFS, Pnoi_dBFS);
      fprintf('#Aircraft Detecting Counting = %d\n', conting);
      disp('======================================');

      if LOG_CSV
        [fid,msg] = fopen(CSV_FILE,'a');
        if fid ~= -1
          fprintf(fid,'%s,%s,%d,%d,%.2f,%.2f,%.2f,%.6f,%.6f,%.1f,%s\n', ...
            ts, ICAO_hex, 17, TC, SNR_iq_dB, Psig_dBFS, Pnoi_dBFS, ...
            lat_deg, lon_deg, range_km, char(src));
          fclose(fid);
        else
          warning('CSV append failed: %s', msg);
        end
      end
    end

  end % st
end % while

%% ---------- Helper: Mode-S/ADS-B CRC-24 ----------
function isValid = checkCRC(bits)
  % Mode S / ADS-B CRC-24 over the 88-bit message field with 24-bit parity.
  % Generator polynomial: 0x1FFF409 (binary 1111111111111010000001001).
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

%% ---------- CPR Helpers: Local CPR robust (no speed check) ----------
% Local CPR using a single frame and a known reference:
% - Try k_ref±1 latitude zones and k_ref±1 longitude zones.
% - Pick the candidate closest to the reference.
% - Reject if range > MAX_RANGE_KM.
function [ok, lat, lon] = cprDecodeLocalRobustSimple(lat_ref, lon_ref, Y, X, F, MAX_RANGE_KM)
  ok = false; lat = NaN; lon = NaN;

  % Latitude grid size depends on frame parity F (even/odd).
  if F == 0, dLat = 360/60; else, dLat = 360/59; end

  % Candidate latitudes around the reference zone index.
  k0 = floor(lat_ref/dLat);
  klist = [k0-1, k0, k0+1];
  latCand = nan(1,3);
  for i=1:3
    lat_i = dLat*(klist(i) + Y);
    if lat_i >= 270, lat_i = lat_i - 360; end     % CPR latitude wrap
    latCand(i) = lat_i;
  end

  % For each latitude candidate, compute NL and the corresponding longitude grid.
  bestDist = Inf; bestLat = NaN; bestLon = NaN;
  for i=1:3
    lat_i = latCand(i);
    NLv = cprNL(lat_i);
    if NLv <= 0, continue; end

    % Longitude grid size depends on NL and F.
    if F==0
      dLon = 360/NLv;
    else
      if NLv <= 1, continue; end
      dLon = 360/(NLv-1);
    end

    % Candidate longitudes around the reference zone index.
    j0 = floor(lon_ref/dLon);
    jlist = [j0-1, j0, j0+1];

    for j=1:3
      lon_j = dLon*(jlist(j) + X);
      if lon_j > 180,  lon_j = lon_j - 360; end
      if lon_j < -180, lon_j = lon_j + 360; end

      dkm = haversine_km(lat_ref, lon_ref, lat_i, lon_j);
      if dkm < bestDist
        bestDist = dkm; bestLat = lat_i; bestLon = lon_j;
      end
    end
  end

  % Guard: must be within a realistic range from the reference.
  if ~isfinite(bestDist) || bestDist > MAX_RANGE_KM
    return;
  end

  ok  = true;
  lat = bestLat;
  lon = bestLon;
end

% DO-260B Annex table: number of longitude zones as a function of latitude.
function NL = cprNL(lat)
  a = abs(lat);
  if     a < 10.47047130, NL = 59;
  elseif a < 14.82817437, NL = 58;
  elseif a < 18.18626357, NL = 57;
  elseif a < 21.02939493, NL = 56;
  elseif a < 23.54504487, NL = 55;
  elseif a < 25.82924707, NL = 54;
  elseif a < 27.93898710, NL = 53;
  elseif a < 29.91135686, NL = 52;
  elseif a < 31.77209708, NL = 51;
  elseif a < 33.53993436, NL = 50;
  elseif a < 35.22899598, NL = 49;
  elseif a < 36.85025108, NL = 48;
  elseif a < 38.41241892, NL = 47;
  elseif a < 39.92256684, NL = 46;
  elseif a < 41.38651832, NL = 45;
  elseif a < 42.80914012, NL = 44;
  elseif a < 44.19454951, NL = 43;
  elseif a < 45.54626723, NL = 42;
  elseif a < 46.86733252, NL = 41;
  elseif a < 48.16039128, NL = 40;
  elseif a < 49.42776439, NL = 39;
  elseif a < 50.67150166, NL = 38;
  elseif a < 51.89342469, NL = 37;
  elseif a < 53.09516153, NL = 36;
  elseif a < 54.27817472, NL = 35;
  elseif a < 55.44378444, NL = 34;
  elseif a < 56.59318756, NL = 33;
  elseif a < 57.72747354, NL = 32;
  elseif a < 58.84763776, NL = 31;
  elseif a < 59.95459277, NL = 30;
  elseif a < 61.04917774, NL = 29;
  elseif a < 62.13216659, NL = 28;
  elseif a < 63.20427479, NL = 27;
  elseif a < 64.26616523, NL = 26;
  elseif a < 65.31845310, NL = 25;
  elseif a < 66.36171008, NL = 24;
  elseif a < 67.39646774, NL = 23;
  elseif a < 68.42322022, NL = 22;
  elseif a < 69.44242631, NL = 21;
  elseif a < 70.45451075, NL = 20;
  elseif a < 71.45986473, NL = 19;
  elseif a < 72.45884545, NL = 18;
  elseif a < 73.45177442, NL = 17;
  elseif a < 74.43893416, NL = 16;
  elseif a < 75.42056257, NL = 15;
  elseif a < 76.39684391, NL = 14;
  elseif a < 77.36789461, NL = 13;
  elseif a < 78.33374083, NL = 12;
  elseif a < 79.29428225, NL = 11;
  elseif a < 80.24923213, NL = 10;
  elseif a < 81.19801349, NL = 9;
  elseif a < 82.13956981, NL = 8;
  elseif a < 83.07199445, NL = 7;
  elseif a < 83.99173563, NL = 6;
  elseif a < 84.89166191, NL = 5;
  elseif a < 85.75541621, NL = 4;
  elseif a < 86.53536998, NL = 3;
  elseif a < 87.00000000, NL = 2;
  else,                       NL = 1;
  end
end

%% ---------- Geo Helper: Great-circle distance (Haversine) ----------
function d_km = haversine_km(lat1, lon1, lat2, lon2)
  rlat1 = deg2rad(lat1); rlon1 = deg2rad(lon1);
  rlat2 = deg2rad(lat2); rlon2 = deg2rad(lon2);
  dlat  = rlat2 - rlat1;
  dlon  = rlon2 - rlon1;
  a     = sin(dlat/2).^2 + cos(rlat1).*cos(rlat2).*sin(dlon/2).^2;
  c     = 2 * atan2(sqrt(a), sqrt(max(0,1-a)));
  R_earth_km = 6371.0;     % mean Earth radius
  d_km  = R_earth_km * c;
end
