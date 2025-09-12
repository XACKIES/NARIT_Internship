%% ADS-B Receiver PlutoSDR ( Model 1 : Direct Xcorr + Adaptive Threshold + Proper Alignment )
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s)
%     -> magnitude (abs)
%     -> matched filter (preamble template, direct xcorr via conv)
%     -> adaptive threshold (median + k*MAD)
%     -> peak pick (candidates above threshold)
%     -> slice after preamble (ALIGN: +L/2 from correlation peak)
%     -> PPM demodulation (bit0/bit1 templates)
%     -> CRC-24 parity check
%
% References:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS).
%   [2] ICAO Annex 10, Vol. IV.
%   [3] Mode S / ADS-B CRC-24 polynomial: 0x1FFF409.
%   [4] Oppenheim & Schafer, Discrete-Time Signal Processing, 3rd ed., Pearson, 2009. (Matched filter, FIR group delay)
%   [5] Huber, Robust Statistics, Wiley, 1981. (median/MAD thresholding)
%
clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT   = false;     % Enable/disable debug plots
DEBUG_EVERY  = 10;        % Plot every N frames
iterCount    = 0;

% Adaptive threshold parameter (k in median + k*MAD)
th_k         = 8.0;      % Typical 6–10. Higher = stricter. (You used 15 before)

%% ---------- Radio / Buffer ----------
fc        = 1090e6;       % ADS-B carrier frequency
sampRate  = 10e6;         % Pluto sampling rate (10 MS/s)
frameLen  = 65536;        % Samples per frame

% PlutoSDR init (Fast Attack AGC)
rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');

%% ---------- ADS-B PHY Parameters ----------
% Oversampling factor relative to 1 Mb/s (1 µs/bit)
os = round(sampRate/1e6);      % ~10 samples per µs @ 10 MS/s

% 16-chip preamble (8 µs); chip = 0.5 µs at 2 MHz chip rate
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];

% --- IMPORTANT: build preamble at the physical chip rate (os/2 samples per 0.5 µs) ---
SamplesPerChip = os/2;         % 1 chip = 0.5 µs → os/2 samples
preamble       = repelem(SyncSequence, SamplesPerChip);  % template used in conv
mf             = conj(flipud(preamble(:)));              % matched-filter kernel

% Bit templates for PPM (1.0 µs per bit = os samples)
half      = floor(os/2);
bit0      = [ones(1,half) zeros(1,os-half)];
bit1      = [zeros(1,half) ones(1,os-half)];

% Pack parameters
ADS_B_Parameter.SamplesPerChip     = SamplesPerChip;           % os/2
ADS_B_Parameter.SamplesPerBit      = os;                       % os samples per 1 µs bit
ADS_B_Parameter.SyncSequence       = SyncSequence;
ADS_B_Parameter.SyncSequenceLength = numel(SyncSequence);
ADS_B_Parameter.PreambleLength     = numel(preamble);          % L = length of template used in conv
ADS_B_Parameter.LongPacketBits     = 112;                      % DF17 total bits (88 payload + 24 parity)

%% ---------- Main Loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Acquire a frame
  rxSig = rx();
  if isempty(rxSig), continue; end

  % 2) Matched filter on magnitude (equivalent to direct xcorr with preamble)
  c = abs(conv(abs(rxSig), mf, 'same'));   % corrOut

  % 3) Adaptive threshold (median + k*MAD)
  med = median(c);
  mad = median(abs(c - med)) + eps;        % robust dispersion (+eps for safety)
  th  = med + th_k * mad;

  % 4) Peak candidates above threshold (simple as your original style)
  locs = find(c > th);
  if isempty(locs), continue; end

  % 5) Candidate extraction and decoding
  for st = locs(:)'   % st is the correlation-peak *near middle of preamble*
    % ALIGNMENT FIX:
    % conv(...,'same') with length L template has group delay ~ L/2,
    % so correlation peak st sits near the *center* of preamble.
    % True payload starts right after preamble end = st + L/2.
    L          = ADS_B_Parameter.PreambleLength;     % template length used in conv
    st_start   = st + floor(L/2);                    % align to end of preamble (start of payload)

    % Guard for buffer end
    msgSamples = ADS_B_Parameter.LongPacketBits * ADS_B_Parameter.SamplesPerBit;
    if st_start + msgSamples - 1 > length(rxSig), continue; end

    % 6) Slice payload region starting at end-of-preamble (properly aligned)
    region = abs(rxSig(st_start : st_start + msgSamples - 1));

    % 7) Vectorized PPM demodulation (bit0/bit1 templates)
    segs   = reshape(region, ADS_B_Parameter.SamplesPerBit, ADS_B_Parameter.LongPacketBits);
    score0 = sum(segs .* bit0.', 1);
    score1 = sum(segs .* bit1.', 1);
    bits   = double(score1 > score0);

    % 8) DF decode and CRC check (Mode S / ADS-B CRC-24)
    DF = bin2dec(num2str(bits(1:5)));
    if DF == 17 && checkCRC(bits)
      CA   = bin2dec(num2str(bits(6:8)));
      ICAO = bits(9:32);
      DATA = bits(33:88);

      disp('========= VALID ADS-B (DF=17) =========');
      fprintf('CA   = %d\n', CA);
      fprintf('ICAO = %06X\n', bin2dec(num2str(ICAO)));
      fprintf('DATA = %s\n', num2str(DATA));
      disp('======================================');

      % Optional debug plot
      if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY) == 0
        try
          figure(1); clf;
          t = (0:numel(c)-1)/sampRate*1e3; % ms
          plot(t, c, 'b-'); grid on; hold on;
          yline(th, '--r', 'th = median + k·MAD');
          % mark st (middle-of-preamble) and st_start (end-of-preamble)
          plot(t(st),       c(st),       'gx', 'MarkerSize', 8, 'LineWidth', 1.5);
          plot(t(st_start), c(st_start), 'mx', 'MarkerSize', 8, 'LineWidth', 1.5);
          legend('corrOut','threshold','st (mid-preamble)','st+L/2 (payload start)','Location','best');
          title('Matched-filter correlation (Direct Xcorr) + Adaptive Threshold + Alignment');
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
