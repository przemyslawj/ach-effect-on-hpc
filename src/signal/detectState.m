function [stateName] = detectState(emgSignal)

dx = diff(emgSignal);



stateName = 'asleep';
end
