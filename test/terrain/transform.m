function b = transform(a, icase)
% =========================================================
%              Generate transformation function
% =========================================================
rand('state', icase);

b = (a - min(min(a))) / (max(max(a)) - min(min(a)));

njump = 4;     % large jumps
pos = rand(njump,1);
mag = rand(njump,1);
shp = rand(njump,1);

b = zeros(size(a));
x = [0:0.01:1];
y = zeros(size(x));
for i = 1:njump
   e = exp((a - pos(i)) * 25 * shp(i));
   b = b + mag(i) * (e - 1./e) ./ (e + 1./e);
   e = exp((x - pos(i)) * 25 * shp(i));
   y = y + mag(i) * (e - 1./e) ./ (e + 1./e);
end
b = (b - y(1)) / (y(end) - y(1));
y = (y - y(1)) / (y(end) - y(1));

nsmall = 10;    % small jumps
pos = rand(nsmall,1);
mag = rand(nsmall,1) / 10;
shp = rand(nsmall,1) * 25;
for i = 1:nsmall
   [m,j] = min(abs(pos(i) - y));
   pos(i) = x(j);
end

for i = 1:nsmall
   e = exp((a - pos(i)) * 25 * shp(i));
   b = b + mag(i) * (e - 1./e) ./ (e + 1./e);
   e = exp((x - pos(i)) * 25 * shp(i));
   y = y + mag(i) * (e - 1./e) ./ (e + 1./e);
end
a = (b - y(1)) / (y(end) - y(1));
y = (y - y(1)) / (y(end) - y(1));

figure; plot(x,y)

b = (a - min(min(a))) / (max(max(a)) - min(min(a)));

