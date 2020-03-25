function y = MovingAverageFilter(b, x)
% Moving average filter calculated with kernel b
    if mod(length(b),2) ~= 1
        error('kernel length should be odd');
    end

    shift = (length(b)-1)/2;

    [y0 z] = filter(b,1,x);

    y = [y0(shift+1:end,:) ; z(1:shift,:)];
end