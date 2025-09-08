%% ADS-B Receiver PlutoSDR ( Model 1 : Direct Xcorr )
% Description:
%   Real-time ADS-B (Mode S, DF=17) receiver using PlutoSDR.
%   Processing chain:
%     PlutoSDR RX (10 MS/s)
%     -> magnitude-squared (energy via abs)
%     -> matched filter (preamble template, direct xcorr)
%     -> peak pick (simple threshold)
%     -> slice after preamble
%     -> PPM demodulation (early-late method via bit0/bit1 templates)
%     -> CRC-24 parity check
%
% References:
%   [1] RTCA DO-260B/DO-260C (1090ES MOPS).
%   [2] ICAO Annex 10, Vol. IV.
%   [3] Mode S / ADS-B CRC-24 polynomial: 0x1FFF409.
%   [4] A. V. Oppenheim and R. W. Schafer, Discrete-Time Signal Processing,
%       3rd ed., Pearson, 2009.  % (Linear-phase FIR / correlation fundamentals)
%
clear; clc;

%% ---------- User Controls ----------
DEBUG_PLOT   = false;     % Enable/disable debug plots (kept for consistency)
DEBUG_EVERY  = 10;        % Plot every N frames
iterCount    = 0;

%% ---------- Radio / Buffer ----------
fc        = 1090e6;       % ADS-B carrier frequency
sampRate  = 10e6;         % Pluto sampling rate
frameLen  = 65536;        % Samples per frame
threshold = 0.01;          % Correlation threshold (relative to peak)
    
% PlutoSDR init (Fast Attack AGC as in the framework)
rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening 1090 MHz ...');


%% ---------- ADS-B PHY Parameters ----------
% Oversampling factor (relative to 1 MHz symbol/bit rate)
os = round(sampRate/1e6);

% 16-chip preamble (8 µs); chip = 0.5 µs at 2 MHz chip rate
SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];

% Matched-filter (direct xcorr) template at sample rate
preamble  = repelem(SyncSequence, os);
mf        = conj(flipud(preamble(:)));

% Bit templates for PPM (1.0 µs per bit = os samples)
half      = floor(os/2);
bit0      = [ones(1,half) zeros(1,os-half)];
bit1      = [zeros(1,half) ones(1,os-half)];

% Pack parameters for readability (consistent with Model 2 style)
ADS_B_Parameter.SamplesPerChip     = os/2;            % os samples/bit -> os/2 per 0.5 µs chip
ADS_B_Parameter.SamplesPerBit      = os;              % samples per 1 µs bit
ADS_B_Parameter.SyncSequence       = SyncSequence;
ADS_B_Parameter.SyncSequenceLength = numel(SyncSequence);
ADS_B_Parameter.PreambleLength     = numel(preamble); % samples
ADS_B_Parameter.LongPacketBits     = 112;             % DF17 payload bits

%% ---------- Main Loop ----------
while true
  iterCount = iterCount + 1;

  % 1) Acquire a frame
  rxSig = rx();
  if isempty(rxSig), continue; end

  % 2) Direct cross-correlation on magnitude (as originally designed)
  corrOut = abs(conv(abs(rxSig), mf, 'same'));

  % 3) Peak picking with relative threshold
  peakVal = max(corrOut);
  locs    = find(corrOut > threshold * peakVal);  % simple threshold detect

  % 4) Candidate extraction and decoding
  for st = locs(:)'
    msgSamples = ADS_B_Parameter.LongPacketBits * ADS_B_Parameter.SamplesPerBit;
    % Slice region immediately after the preamble peak index "st"
    if st + msgSamples - 1 > length(rxSig), continue; end

    region = abs(rxSig(st : st + msgSamples - 1));

    % 5) Vectorized PPM demodulation using bit templates (no loop)
    segs   = reshape(region, ADS_B_Parameter.SamplesPerBit, ADS_B_Parameter.LongPacketBits);
    score0 = sum(segs .* bit0.', 1);
    score1 = sum(segs .* bit1.', 1);
    bits   = double(score1 > score0);

    % 6) DF decode and CRC check
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

      % Optional debug plot (kept style-consistent; does not alter flow)
      if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY) == 0
        try
          figure(1); clf;
          t = (0:numel(corrOut)-1)/sampRate*1e3; % ms
          plot(t, corrOut); grid on;
          title('Matched-filter correlation (Direct Xcorr)');
          xlabel('Time (ms)'); ylabel('Correlation');
          hold on; plot(t(st), corrOut(st), 'rx');
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
