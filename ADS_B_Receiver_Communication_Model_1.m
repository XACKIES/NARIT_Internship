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
%   [4] Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd ed., Pearson, 2009.
%       (Matched filter; FIR linear-phase & group delay).
%   [5] Huber, Robust Statistics, Wiley, 1981. (Median/MAD thresholding).
%   [6] Lyons, Understanding Digital Signal Processing, 3rd ed., Prentice Hall, 2010.
%       (Streaming FIR, passband gain, practical group-delay handling).
%   [7] MathWorks: filter (Initial & final conditions z_i/z_f) – streaming FIR usage.
%
% Notes:
%   - Requires: Communications Toolbox Support Package for ADALM-PLUTO Radio.
%   - AGC Fast Attack recommended at 10 MS/s for bursty ADS-B.
%   - All time alignment is performed on the FIR-filtered stream to keep st+L/2 valid.

clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT   = false;      % Enable/disable debug plots
DEBUG_EVERY  = 10;         % Plot every N frames (when DEBUG_PLOT = true)
iterCount    = 0;

% Adaptive threshold parameter (k in median + k*MAD)
th_k         = 5.0;        % Typical 6–10. Higher = stricter.

%% ---------- Radio / Buffer ----------
fc        = 1090e6;        % ADS-B carrier frequency
sampRate  = 10e6;          % Pluto sampling rate (10 MS/s)
frameLen  = 65536;         % Samples per frame

% PlutoSDR init (Fast Attack AGC)
rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');
conting = 0;
%% ---------- FIR front-end (pre-buffer) ----------
% Linear-phase FIR (real, symmetric). Sum(b_fir) ~ 1 → ~unity DC/passband gain.
b_fir = [0  -0.001951 -0.001727 0.002622 0.006403 0 -0.013408 -0.011526 0.015748 ...
         0.034483 0 -0.064286 -0.056995 0.089999 0.300257 0.400761 0.300257 ...
         0.089999 -0.056995 -0.064286 0 0.034483 0.015748 -0.011526 -0.013408 ...
         0 0.006403 0.002622 -0.001727 -0.001951 0];

M_fir  = numel(b_fir);
D_fir  = (M_fir-1)/2;                  % group delay in samples (31 taps ⇒ 15 samples ≈ 1.5 µs @10 MS/s)

% Streaming filter state (z_i) for FIR
zi_fir = zeros(M_fir-1,1);

%% ---------- ADS-B PHY Parameters ----------
% Oversampling factor relative to 1 Mb/s (1 µs/bit)
os = round(sampRate/1e6);              % ~10 samples per µs @ 10 MS/s

% 16-chip preamble (8 µs); chip = 0.5 µs at 2 MHz chip rate
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];

% Build preamble at the physical chip rate (os/2 samples per 0.5 µs)
SamplesPerChip = max(1, round(os/2));                    % integer replication count
preamble       = repelem(SyncSequence, SamplesPerChip);  % template used in conv
mf             = flipud(preamble(:));                    % matched-filter kernel (real)

% Bit templates for PPM (1.0 µs per bit = os samples)
half      = floor(os/2);
bit0      = [ones(1,half) zeros(1,os-half)];
bit1      = [zeros(1,half) ones(1,os-half)];

% Pack parameters (for docs/telemetry)
ADS_B_Parameter = struct();
ADS_B_Parameter.SamplesPerChip     = SamplesPerChip;           % ≈ os/2
ADS_B_Parameter.SamplesPerBit      = os;                       % os samples per 1 µs bit
ADS_B_Parameter.SyncSequence       = SyncSequence;
ADS_B_Parameter.SyncSequenceLength = numel(SyncSequence);
ADS_B_Parameter.PreambleLength     = numel(preamble);          % L = length of template used in conv
ADS_B_Parameter.LongPacketBits     = 112;                      % DF17 total bits (88 payload + 24 parity)
ADS_B_Parameter.FIRCoeffs          = b_fir;
ADS_B_Parameter.FIRLength          = M_fir;
ADS_B_Parameter.FIRGroupDelay      = D_fir;                    % samples

%% ---------- Main Loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Acquire a frame
  rxSig = rx();
  if isempty(rxSig), continue; end

  % 2) Front-end FIR on complex IQ (streaming with state)
  [y_fir, zi_fir] = filter(b_fir, 1, rxSig, zi_fir);   % complex-in/complex-out

  % 3) Matched filter on magnitude (equivalent to direct xcorr with preamble)
  %    Use filtered stream for the entire downstream chain.
  c = abs(conv(abs(y_fir), mf, 'same'));   % corrOut

  % 4) Adaptive threshold (median + k*MAD)
  med = median(c);
  mad = median(abs(c - med)) + eps;        % robust dispersion (+eps for safety)
  th  = med + th_k * mad;

  % 5) Peak candidates above threshold (simple selection)
  locs = find(c > th);
  if isempty(locs), continue; end

  % 6) Candidate extraction and decoding
  for st = locs(:)'   % st is the correlation-peak *near middle of preamble*
    % ALIGNMENT:
    % conv(...,'same') with length L template has group delay ~ L/2,
    % so correlation peak st sits near the *center* of preamble.
    % True payload starts right after preamble end = st + L/2.
    L          = ADS_B_Parameter.PreambleLength;      % template length used in conv
    st_start   = st + floor(L/2);                     % align to end of preamble (start of payload)

    % Guard for buffer end
    msgSamples = ADS_B_Parameter.LongPacketBits * ADS_B_Parameter.SamplesPerBit;
    if st_start < 1 || st_start + msgSamples - 1 > length(y_fir), continue; end

    % 7) Slice payload region starting at end-of-preamble (on filtered magnitude)
    region = abs(y_fir(st_start : st_start + msgSamples - 1));

    % 8) Vectorized PPM demodulation (bit0/bit1 templates)
    segs   = reshape(region, ADS_B_Parameter.SamplesPerBit, ADS_B_Parameter.LongPacketBits);
    score0 = sum(segs .* bit0.', 1);
    score1 = sum(segs .* bit1.', 1);
    bits   = double(score1 > score0);

    % 9) DF decode and CRC check (Mode S / ADS-B CRC-24)
    DF = bin2dec(char(bits(1:5) + '0'));   % robust: no spaces
    if DF == 17 && checkCRC(bits)
      CA   = bin2dec(char(bits(6:8)    + '0'));
      ICAO = bits(9:32);
      DATA = bits(33:88);
      conting = conting+1;
      disp('========= VALID ADS-B (DF=17) =========');   
      fprintf('CA   = %d\n', CA);
      fprintf('ICAO = %06X\n', bin2dec(char(ICAO + '0')));
      fprintf('DATA = %s\n', char(DATA + '0'));
      fprintf('#Aircraft Detecting Counting = %d\n', (conting));
      disp('======================================');
      % Optional debug plot
      if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY) == 0
        try
          figure(1); clf;
          t = (0:numel(c)-1)/sampRate*1e3; % ms
          plot(t, c, '-'); grid on; hold on;
          yline(th, '--', 'th = median + k·MAD');
          % mark st (middle-of-preamble) and st_start (end-of-preamble)
          plot(t(st),       c(st),       'gx', 'MarkerSize', 8, 'LineWidth', 1.5);
          plot(t(st_start), c(st_start), 'mx', 'MarkerSize', 8, 'LineWidth', 1.5);
          legend('corrOut','threshold','st (mid-preamble)','st+L/2 (payload start)','Location','best');
          title('Matched-filter correlation (Direct Xcorr) + FIR + Adaptive Threshold + Alignment');
          xlabel('Time (ms)'); ylabel('Correlation magnitude');
          drawnow limitrate;
        catch
        end
      end
    end
  end
end

%% ---------- Helper: Mode-S/ADS-B CRC-24 parity ----------
function isValid = checkCRC(bits)
  % Mode S / ADS-B (DO-260) CRC-24 over 88-bit payload with 24 parity bits.
  % Generator polynomial (binary): 1111111111111010000001001  = 0x1FFF409.
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
