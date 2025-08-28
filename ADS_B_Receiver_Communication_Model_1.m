%% ADS-B Receiver PlutoSDR
clear; 

%% Parameters
fc       = 1090e6;       % ADS-B frequency
sampRate = 10e6;          % Pluto sampling rate (2 MHz พอแล้ว)
frameLen = 65536;        % Samples per frame (ใหญ่ = เร็วขึ้น)
thresh   = 0.1;         % Correlation threshold

%% PlutoSDR init
rx = sdrrx('Pluto', ...
    'CenterFrequency', fc, ...
    'BasebandSampleRate', sampRate, ...
    'SamplesPerFrame', frameLen, ...
    'GainSource','Manual', ...
    'Gain', 50, ...
    'OutputDataType','double');

disp('✅ PlutoSDR ready... Listening...');

%% Oversampling factor
os = round(sampRate/1e6);

%% Preamble sequence (8 µs = 16 chips)
preambleBits = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];
preamble = repelem(preambleBits, os);
mf = conj(flipud(preamble(:)));

%% Templates for bit0/bit1 (PPM)
half = floor(os/2);
bit0 = [ones(1,half) zeros(1,os-half)];
bit1 = [zeros(1,half) ones(1,os-half)];

%% Main Loop
while true
    rxSig = rx(); 
    if isempty(rxSig), continue; end
    
    % === Step 1: Correlation ===
    corrOut = abs(conv(abs(rxSig), mf, 'same'));
    peakVal = max(corrOut);
    locs = find(corrOut > thresh*peakVal);  % simple threshold detect
    
    % === Step 2: Extract candidate packets ===
    for st = locs(:)'
        msgSamples = 112*os;
        if st+msgSamples-1 > length(rxSig), continue; end
        
        region = abs(rxSig(st:st+msgSamples-1));
        
        % Vectorized reshape -> no loop
        segs = reshape(region, os, 112);
        score0 = sum(segs .* bit0.', 1);
        score1 = sum(segs .* bit1.', 1);
        bits = double(score1 > score0);
        
        % === Step 3: Decode ===
        DF = bin2dec(num2str(bits(1:5)));
        if DF==17 && checkCRC(bits)
            CA   = bin2dec(num2str(bits(6:8)));
            ICAO = bits(9:32);
            DATA = bits(33:88);
            
            disp('========= VALID ADS-B (DF=17) =========');
            fprintf('CA   = %d\n', CA);
            fprintf('ICAO = %06X\n', bin2dec(num2str(ICAO)));
            fprintf('DATA = %s\n', num2str(DATA));
            disp('======================================');
        end
    end
end

%% CRC Function (24-bit)
function isValid = checkCRC(bits)
    poly = [1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1];  
    msg   = bits(1:88);
    crcRx = bits(89:112);
    dividend = [msg zeros(1,24)];
    for i=1:length(msg)
        if dividend(i)==1
            dividend(i:i+24) = bitxor(dividend(i:i+24), poly);
        end
    end
    crcCalc = dividend(end-23:end);
    isValid = isequal(crcCalc, crcRx);
end
