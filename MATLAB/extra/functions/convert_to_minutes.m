function dt = convert_to_minutes(duration)
    [hours, minutes, seconds] = hms(duration);
    dt = hours * 60 + minutes + seconds / 60;
end
