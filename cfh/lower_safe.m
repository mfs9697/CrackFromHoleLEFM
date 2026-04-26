function s = lower_safe(x)
    if isstring(x)
        x = char(x);
    end
    if ~ischar(x)
        error('build_appended_hole_loop:BadStringInput', ...
            'Expected char or string input.');
    end
    s = x;
    mask = (s >= 'A') & (s <= 'Z');
    s(mask) = char(double(s(mask)) + ('a' - 'A'));
end