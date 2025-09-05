%% ADS-B Receiver PlutoSDR (Final, Clean Version)
clear; clc;

DEBUG_PLOT   = false;
DEBUG_EVERY  = 10;
iterCount    = 0;

fc       = 1090e6;
sampRate = 10e6;
frameLen = 65536;
bufferLen= 10*frameLen;

rx = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'BasebandSampleRate', sampRate, ...
  'SamplesPerFrame', frameLen, ...
  'GainSource','AGC Fast Attack', ...
  'OutputDataType','double');

disp('✅ PlutoSDR ready (AGC Fast Attack) ... Listening 1090 MHz ...');

spc_chip  = round(sampRate/2e6);
spb_bit   = 2*spc_chip;

SyncSequence = [1 0 1 0 0 0 0 1 0 1 0 0 0 0 0 0];

adsbParam.SamplesPerChip          = spc_chip;
adsbParam.SyncSequence            = SyncSequence;
adsbParam.SyncSequenceLength      = length(SyncSequence);
adsbParam.SyncSequenceHighIndices = find(SyncSequence==1);
adsbParam.SyncSequenceLowIndices  = find(SyncSequence==0);
adsbParam.SyncSequenceNumHighValues = numel(adsbParam.SyncSequenceHighIndices);
adsbParam.SyncSequenceNumLowValues  = numel(adsbParam.SyncSequenceLowIndices);
adsbParam.PreambleLength          = adsbParam.SyncSequenceLength * spc_chip;
adsbParam.LongPacketLength        = 112 * spb_bit;
adsbParam.MaxPacketLength         = adsbParam.PreambleLength + adsbParam.LongPacketLength;
adsbParam.MaxNumPacketsInFrame    = 64;
adsbParam.SyncDownsampleFactor    = 1;

preamble_bip = 2*SyncSequence - 1;
mf = flipud(repelem(preamble_bip(:), spc_chip));

xBuff = zeros(bufferLen,1);

while true
  iterCount = iterCount + 1;
  newFrame = rx(); if isempty(newFrame), continue; end
  xBuff = [xBuff(numel(newFrame)+1:end); newFrame];
  energySig = abs(xBuff).^2;
  xFilt = conv(energySig, mf, 'same');
  [packetSamples, packetCnt, syncTimeVec] = packetSearch(abs(xFilt), xBuff, energySig, adsbParam);

  if DEBUG_PLOT && mod(iterCount, DEBUG_EVERY)==0
    try
      figure(1); clf;
      t = (1:numel(energySig))/sampRate*1e3;
      subplot(2,1,1);
      plot(t, energySig); grid on;
      hold on;
      if ~isempty(syncTimeVec)
        plot(t(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), ...
             energySig(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), 'ro');
      end
      subplot(2,1,2);
      plot(t, xFilt); grid on;
      hold on;
      if ~isempty(syncTimeVec)
        plot(t(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), ...
             xFilt(syncTimeVec(syncTimeVec>0 & syncTimeVec<=numel(t))), 'rx');
      end
      drawnow limitrate;
    catch
    end
  end

  for k = 1:packetCnt
    region = abs(packetSamples(:,k));
    nBits = floor(numel(region)/spb_bit); if nBits < 112, continue; end
    segs  = reshape(region(1:nBits*spb_bit), spb_bit, nBits);
    early = sum(segs(1:spc_chip, :), 1);
    late  = sum(segs(spc_chip+1:end, :), 1);
    bits  = double(early > late);
    bits  = bits(1:112);

    if checkCRC(bits)
      DF = bin2dec(num2str(bits(1:5)));
      if DF==17
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
end

function [packetSamples, packetCnt, syncTimeVec] = packetSearch(xFilt, xBuff, energySig, adsbParam)
  spc = adsbParam.SamplesPerChip;
  syncLen    = adsbParam.SyncSequenceLength;
  syncSigLen = syncLen*spc;
  xLen       = length(xBuff);

  subFrameLen     = adsbParam.MaxPacketLength;
  subFrameDownLen = subFrameLen / adsbParam.SyncDownsampleFactor;
  numSubFrames    = floor(xLen / subFrameLen);

  packetSamples = zeros(adsbParam.LongPacketLength, adsbParam.MaxNumPacketsInFrame, 'like', xBuff);
  syncTimeVec   = zeros(adsbParam.MaxNumPacketsInFrame,1);
  packetCnt     = 0;

  for p = 0:(numSubFrames-2)
    idx = double(p)*subFrameDownLen + (1:subFrameDownLen);
    [~, tmp] = max(xFilt(idx));
    syncIdx  = tmp;
    syncTime = round(syncIdx*adsbParam.SyncDownsampleFactor - syncSigLen + p*subFrameLen);

    if (syncTime <= 0) || (syncTime + adsbParam.MaxPacketLength - 1 > xLen)
      continue;
    end

    rxSyncEnergy = energySig(syncTime + (0:syncSigLen-1));
    rxSyncSeq    = sum(reshape(rxSyncEnergy, spc, syncLen), 1);
    hi = adsbParam.SyncSequenceHighIndices;
    lo = adsbParam.SyncSequenceLowIndices;
    th = (mean(rxSyncSeq(hi)) + mean(rxSyncSeq(lo)))/2;

    if all(xor((rxSyncSeq < th), adsbParam.SyncSequence))
      packetCnt = packetCnt + 1;
      if packetCnt <= adsbParam.MaxNumPacketsInFrame
        dataIdx = int32(adsbParam.PreambleLength + (0:adsbParam.LongPacketLength-1));
        packetSamples(:,packetCnt) = xBuff(syncTime + dataIdx, 1);
        syncTimeVec(packetCnt)     = syncTime;
      end
    end
  end
end

function isValid = checkCRC(bits)
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
