% =========================================================
%              Generate random image
% =========================================================
pertsize = [1.0, 2.0, 4.0, 6.0, 4.0, 2.0, 1.0, 0.5, 0.2];
niter = 8;

% rand('state', 100); itrans = 1; sea = 0.075;   % Mont Olympus
rand('state', 101); itrans = 2; sea = 0.0;     % North - South continent
% rand('state', 102); itrans = 5; sea = 0.20;      % Tharsis - plains

a0 = rand([3,3]);
for iter = 1:niter
   a0(1,:) = 0.0; a0(end,:) = 0.0; a0(:,1) = 0.0; a0(:,end) = 0.0;

   newsize = [size(a0,1)*2-1, size(a0,2)*2-1];

   p = rand(newsize) - 0.5;
   for i = 1:0
      p(2:2:end,:) = (p(1:2:end-1,:) + 2*p(2:2:end,:) + p(3:2:end,:)) / 4;
   end

   a = zeros(newsize);
   a(1:2:end, 1:2:end) = a0;
   a(1:2:end, 2:2:end) = (a(1:2:end, 1:2:end-1) + a(1:2:end, 3:2:end)) / 2.0;
   a(2:2:end, 1:2:end) = (a(1:2:end-1, 1:2:end) + a(3:2:end, 1:2:end)) / 2.0;
   a(2:2:end, 2:2:end) = (a(2:2:end, 1:2:end-1) + a(2:2:end, 3:2:end)) / 2.0;
   a = a + p * 0.5^iter * pertsize(iter);
   a(1:2:end, 1:2:end) = a0;
   a0 = a;
end

a(:,:) = (a(:,:) - min(min(a(:,:)))) / (max(max(a(:,:))) - min(min(a(:,:))));

% =========================================================
%              Generate transformation function
% =========================================================
a = transform(a,itrans);

% =========================================================
%              Draw image
% =========================================================
% a(:,:,1) = a(:,:,1) * 0.3 + 0.1;
% a(:,:,2) = a(:,:,2) * 0.6 + 0.3;
% a(:,:,3) = a(:,:,3) * 0.3 + 0.1;
if sea < a(1,1)
   sea = a(1,1) + 0.01;
end
figure;

% contour plot
a = (a - sea) / (1 - sea);
a(a < 0) = - 0.7;
contour(a,25);
axis square;

% surf plot
% a(a == 0) = 0.2 / (1.2-sea);
% surf(a*2);
% axis equal;
