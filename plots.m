clear;
close all;

% displaying the image P
P_orig = load('P.mat');
P_orig = P_orig.P; 
P_orig = P_orig/10;    % scaling the array to cm for radon transform  
figure;
imshow(P_orig,[]);
title("phantom image P");   

% taking the radon transform of the image 
angles = 0:1:179;
[g, xp] = radon(P_orig,angles);  

% taking inverse radon transform
P_rec = iradon(g,angles);
figure;
imshow(P_rec,[]);
title("recontructed image P without noise"); 

% adding the poisson noise to the projection 
N0 = 10^3;
[m,n] = size(g);
g_free = g(:,1);
N = N0 * exp(-g_free);
r = poissrnd(N);
N = r;
g_noise = -log(0.0001 + N/N0);      % adding a very small term to avoid infinite value
                                    % does not make a change since to small 
%%
% plotting the noise free and noise added projections in image domain
FOV = 51.2;
b = (length(xp)-1)/2;
l = ( -b:1:b )* (FOV/(2*b));
figure;plot(l,g_free)
title("projection without noise");xlabel("( cm )");
figure;plot(l,g_noise)
title("projection with noise added");xlabel("( cm )");

% plotting the noise free and noise added projections in fourier domain
f = ( -b:1:b )* ((2*b)/FOV);
G_free = fftshift(fft(g_free));
figure; plot(f,abs(G_free));
title("noise free projection in fourier domain");xlabel("( cm-1 )");
G_noise = fftshift(fft(g_noise));
figure; plot(f,abs(G_noise));
title("noisy projection in fourier domain");xlabel("( cm-1 )");

%%
% creating a windowed ramp filter 
freqs=linspace(-1, 1, length(xp))' ; 
ramp_filter = abs(freqs);
ramp_filter = repmat(ramp_filter, [1 length(angles)]);

% considering the detectors as having width d 
d = 0.1;   
V = sym(d * freqs);
detector_filter = d * sinc(V) ;
detector_filter = repmat(detector_filter, [1 length(angles)]);
detector_filter = double(detector_filter);

% plotting the filters in fourier domain
figure;plot(f,ramp_filter);
title("ramp filter in fourier domain");xlabel("cm(-1)");
figure; plot(f,detector_filter)
title("detector filter in fourier domain");xlabel("cm(-1)");
