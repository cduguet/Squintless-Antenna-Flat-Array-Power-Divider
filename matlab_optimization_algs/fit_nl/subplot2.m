function subplot2(m, n, p)
% Similar to subplot, but p can exceed the number of subplots per figure
% (rows * columns = m * n). A new figure will then be created automatically.
% Always start with p == 1 or a multiple of rows * columnms or create figure
% before.
%
% Syntax
%   function subplot2(m, n, p)

d = mod(p, m * n);

if d == 0
    d = m * n;
end

if d == 1
    figure;
end

subplot(m, n, d);
end
