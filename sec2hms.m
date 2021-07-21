function hms = sec2hms(t)
    hours = floor(t / 3600);
    t = t - hours * 3600;
    mins = floor(t / 60);
    secs = t - mins * 60;
    hms = sprintf('%02d:%02d:%05.2f\n', hours, mins, secs);
end
