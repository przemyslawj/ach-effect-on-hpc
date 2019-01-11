function [dataArray] = filter50Hz(dataArray, fs)
% Filter 50 Hz with IIR notch filter

wo = 50 / (fs / 2);  bw = wo / 35;
[b, a] = iirnotch(wo,bw);

for i = 1:size(dataArray, 1)
    dataArray(i,:) = filter(b,a, dataArray(i,:));
end

end
