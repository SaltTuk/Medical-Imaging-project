clear;
close all;

% create the attenuation coefficient transform polynomials
p = mu_polyfit;
P_orig = createPhantom(p);

%% part 1: bandwith effect on image noise variance  
bw = [500 1000 1500 2000 2500];
P_rec_noise1 = zeros(5,514,514);
av_var1 = zeros(1,5);
for i = 1:5
    [P_rec_noise1(i,:,:), av_var1(i)] = bandwith_effect(P_orig,bw(i));
end

figure;
subplot(2,1,1);imshow(P_rec_noise1(1,:,:), []); 
title(sprintf('FBP with %.0f cm-1 bandwith rectangular window', bw(1)));
xlabel("average variance = " + av_var1(1));
 
subplot(2,1,2);imshow(P_rec_noise1(5,:,:), []);
title(sprintf('FBP with %.0f cm-1 bandwith rectangular window', bw(5)));
xlabel("average variance = " + av_var1(5));

%% part 2: detector width effect on image noise variance  
d = [2 4 6 8 10];
P_rec_noise2 = zeros(5,514,514);
av_var2 = zeros(1,5);
for i = 1:5
    [P_rec_noise2(i,:,:), av_var2(i)] = detector_effect(P_orig,d(i));
end

figure;
subplot(2,1,1);imshow(P_rec_noise2(1,:,:), []); 
title(sprintf('FBP with %.0f mm detector width', 10*d(1)));
xlabel("average variance = " + av_var2(1));
 
subplot(2,1,2);imshow(P_rec_noise2(5,:,:), []);
title(sprintf('FBP with %.0f mm detector width', 10*d(5)));
xlabel("average variance = " + av_var2(5));

%% part 3: incident photon count effect on image noise variance  
N0 = [2 4 6 8 10]*10^5;
P_rec_noise3 = zeros(5,514,514);
av_var3 = zeros(1,5);
for i = 1:5
    [P_rec_noise3(i,:,:), av_var3(i)] = photon_effect(P_orig,N0(i));
end


figure;
subplot(2,1,1);imshow(P_rec_noise3(1,:,:), []); 
title(sprintf('FBP with %.0f incident photon number', N0(1)));
xlabel("average variance = " + av_var3(1));
 
subplot(2,1,2);imshow(P_rec_noise3(5,:,:), []);
title(sprintf('FBP with %.0f incident photon number', N0(5)));
xlabel("average variance = " + av_var3(5));

%% part 4: number of projection angles effect on image noise variance 
angle_num = [100 187 275 362 450];
P_rec_noise4 = zeros(5,514,514);
av_var4 = zeros(1,5);
for i = 1:5
    [P_rec_noise4(i,:,:), av_var4(i)] = num_angle_effect(P_orig,angle_num(i));
end

figure;
subplot(2,1,1);imshow(P_rec_noise4(1,:,:), []); 
title(sprintf('FBP with %.0f number of projection angles', angle_num(1)));
xlabel("average variance = " + av_var4(1));
 
subplot(2,1,2);imshow(P_rec_noise4(5,:,:), []);
title(sprintf('FBP with %.0f number of projection angles', angle_num(5)));
xlabel("average variance = " + av_var4(5));

%% FUNCTIONS 

% part 1: bandwith effect on image noise variance   
function [P_rec_noise, av_var] = bandwith_effect(P_orig,bandwith)

    % taking the radon transform of the image 
    angles = 0:1:179;
    [g, xp] = radon(P_orig,angles);  
    P_run = zeros(10,514,514);
    % 10 independent measurements
    for sa = 1:10
        % adding the poisson noise to the projections at each angle 
        N0 = poissrnd(10^5);
        [m,n] = size(g);
        gN = zeros(m,n);
        for i = 1:n
            N = N0 * exp(-g(:,i));
            r = poissrnd(N);
            N = r;
            gN(:,i) = -log(0.0001 + N/N0);   % adding a very small term to avoid infinite value
        end                                 % does not make a change since to small 

        % applying rectangular window to the ramp filter 
        % creating a windowed ramp filter 
        freqs=linspace(-1, 1, length(xp))' ; 
        ramp_filter = abs(freqs);

        % applying the bandwith to the signals in the fourier domain
        FOV = 51.2;
        b = (length(xp)-1)/2;
        bw = (2*b^2)/FOV;
        index = round(((bw-bandwith)/bw)*183);
        ramp_filter(1:index) = 0;
        ramp_filter(length(xp)-index:length(xp)) = 0;
        ramp_filter = repmat(ramp_filter, [1 length(angles)]);

        % considering the detectors as having width d 
        d = 0.1;   
        V = sym(d * freqs);
        detector_filter = d * sinc(V) ;
        detector_filter = repmat(detector_filter, [1 length(angles)]);
        detector_filter = double(detector_filter);

        % doing FT domain filtering 
        G = fftshift(fft(gN,[],1),1);
        G_filtered = G .* ramp_filter .* detector_filter;

        % taking inverse transform 
        g_filtered = ifftshift(G_filtered,1);
        g_filtered = real(ifft(g_filtered,[],1));

        % FBP without any filter since we have applied own filter already 
        P_rec_noise = iradon(g_filtered,angles,"linear", 'None');

        % calculating the image noise variance 
        P_run(sa,:,:) = P_rec_noise;
    end
    % calculate the variance among 10 measurements, and the average
    % variance among the elements
    av_var = mean(mean(var(P_run)));
end

% part 2: detector width effect on image noise variance  
function [P_rec_noise, av_var] = detector_effect(P_orig,detector_width)
    
    % taking the radon transform of the image 
    angles = 0:1:179;
    [g, xp] = radon(P_orig,angles);  
    P_run = zeros(10,514,514);
    
    % 10 independent measurements
    for sa = 1:10
        % adding the poisson noise to the projections at each angle 
        N0 = poissrnd(10^5);
        [m,n] = size(g);
        gN = zeros(m,n);
        for i = 1:n
            N = N0 * exp(-g(:,i));
            r = poissrnd(N);
            N = r;
            gN(:,i) = -log(0.0001 + N/N0);   % adding a very small term to avoid infinite value
        end                                 % does not make a change since to small 

        % applying rectangular window to the ramp filter 
        % creating a windowed ramp filter 
        freqs=linspace(-1, 1, length(xp))' ; 
        ramp_filter = abs(freqs);
        ramp_filter = repmat(ramp_filter, [1 length(angles)]);

        % considering the detectors as having width d 
        d = detector_width;   
        V = sym(d * freqs);
        detector_filter = d * sinc(V) ;
        detector_filter = repmat(detector_filter, [1 length(angles)]);
        detector_filter = double(detector_filter);

        % doing FT domain filtering 
        G = fftshift(fft(gN,[],1),1);
        G_filtered = G .* ramp_filter .* detector_filter;  % adding the detector window effect
        g_filtered = ifftshift(G_filtered,1);
        g_filtered = real(ifft(g_filtered,[],1));

        % FBP without any filter since we have applied own filter already 
        P_rec_noise = iradon(g_filtered,angles,"linear", 'None');

        % calculating the image noise variance 
        P_run(sa,:,:) = P_rec_noise;
    end
    % calculate the variance among 10 measurements, and the average
    % variance among the elements
    av_var = mean(mean(var(P_run)));
end

% part 3: incident photon count effect on image noise variance  
function [P_rec_noise, av_var] = photon_effect(P_orig,photon_count)
    % taking the radon transform of the image 
    angles = 0:1:179;
    [g, xp] = radon(P_orig,angles);   
    P_run = zeros(10,514,514);
    
    % 10 independent measurements
    for sa = 1:10
        % adding the poisson noise to the projections at each angle 
        N0 = poissrnd(photon_count);
        [m,n] = size(g);
        gN = zeros(m,n);
        for i = 1:n
            N = N0 * exp(-g(:,i));
            r = poissrnd(N);
            N = r;
            gN(:,i) = -log(0.0001 + N/N0);   % adding a very small term to avoid infinite value
        end                                 % does not make a change since to small 

        % applying rectangular window to the ramp filter 
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

        % doing FT domain filtering 
        G = fftshift(fft(gN,[],1),1);
        G_filtered = G .* ramp_filter .* detector_filter;
        g_filtered = ifftshift(G_filtered,1);
        g_filtered = real(ifft(g_filtered,[],1));

        % FBP without any filter since we have applied own filter already 
        P_rec_noise = iradon(g_filtered,angles,"linear", 'None');

        % calculating the image noise variance 
        P_run(sa,:,:) = P_rec_noise;
    end
    % calculate the variance among 10 measurements, and the average
    % variance among the elements
    av_var = mean(mean(var(P_run)));
end

% part 4: number of projection angles effect on image noise variance   
function [P_rec_noise, av_var] = num_angle_effect(P_orig,num_angle)

    % taking the radon transform of the image 
    d_angle = 180 / num_angle;
    angles = 0:d_angle:179;
    [g, xp] = radon(P_orig,angles); 
    
    % 10 independent measurements
    P_run = zeros(10,514,514);
    for sa = 1:10
        % adding the poisson noise to the projections at each angle 
        N0 = poissrnd(10^5);
        [m,n] = size(g);
        gN = zeros(m,n);
        for i = 1:n
            N = N0 * exp(-g(:,i));
            r = poissrnd(N);
            N = r;
            gN(:,i) = -log(0.0001 + N/N0);   % adding a very small term to avoid infinite value
        end                                 % does not make a change since to small 

        % applying rectangular window to the ramp filter 
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

        % doing FT domain filtering 
        G = fftshift(fft(gN,[],1),1);
        G_filtered = G .* ramp_filter .* detector_filter;
        g_filtered = ifftshift(G_filtered,1);
        g_filtered = real(ifft(g_filtered,[],1));

        % FBP without any filter since we have applied own filter already 
        P_rec_noise = iradon(g_filtered,angles,"linear", 'None');

        % calculating the image noise variance 
        P_run(sa,:,:) = P_rec_noise;
    end
    % calculate the variance among 10 measurements, and the average
    % variance among the elements
    av_var = mean(mean(var(P_run)));
end

% create and display the transformed phantom from chest CT
function [P] = createPhantom(p)
    body_tissues = {'outside phantom','adrenals','skin','fat','spinal cord','blood pool','spine','gas (bowel)','rib cage & sternum','fluid (bowel)','pelvis','bone marrow','long bones','lymph nodes','skeletal muscle','thyroid','lungs','trachea','heart','cartilage','liver','spleen','gall bladder','urine','kidney','feces','pharynx','testes','esophagus','prostate','stomach','rectum','small bowel','diaphragm','colon','bladder','pancreas','lesion'};
    bt_ids = [0,84,64,85,66,86,68,87,69,88,70,89,71,90,72,91,73,92,74,93,75,94,76,95,77,96,78,97,79,98,80,100,81,102,82,103,83,126];
    head_tissues = {'outside phantom','medulla oblongata','skin','fat','brain','blood pool','spinal','bone marrow','skull','pons','spine','trachea','dens of axis','cartilage','jaw bone','uncus (ear bones)','skeletal muscle','sinuses/mouth cavity','lachrymal glands','optic nerve','spinal canal','cerebral falx','hard palate','eye','cerebellum','lens','tongue','cerebral aqueduct','pharynx','teeth','esophagus'};
    ht_ids = [0,84,64,85,65,86,66,89,67,90,68,92,69,93,70,98,72,103,73,105,74,112,75,118,76,120,77,121,78,124,79];

    % choose the phantom slice
    segment = 'BODY';
    slice = 15;

    % read color segmented image file convert
    fileID = fopen(strcat('dataset\Color_slices\',segment,'\',segment,'_',num2str(slice,'%03d'),'_C.DAT.1'));
    P_clr = fread(fileID,[512,512],'uint8=>double').';

    % read CT image file convert
    fileID = fopen(strcat('dataset\CT_slices\',segment,'\',segment,'_',num2str(slice,'%03d'),'_I.DAT.1'));
    P_ct = fread(fileID,[512,512],'int16=>int16');
    P_ct = double(swapbytes(P_ct)).';
    P_ct_mu = P_ct*p(1) + p(2); % apply transform polynomials
    P_ct_mu(P_clr == 0 | P_ct_mu < 0) = 0; % mask the air regions and the negative values

    % find the unique ids in the segmented image
    [Aval,~,indAval] = unique(P_clr);

    % create an array of the names of those ids
    if strcmp(segment,'BODY')
        [~,ia,~] = intersect(bt_ids,Aval);
        tissues = body_tissues(ia);
    else
        [~,ia,~] = intersect(ht_ids,Aval);
        tissues = head_tissues(ia);
    end

    % create a new segmented image with ids replaced with 1,2,3,...
    Avalnew = 1:length(Aval); 
    Anew = Avalnew(indAval);
    Anew = reshape(Anew, size(P_clr));

    % plot the image with a distinct color for each region
    cmap = turbo(length(Aval));
    cmap(1,:) = [0,0,0];
    figure(1)
    subplot(1,3,2)
    image(Anew);
    colormap(cmap);
    hold on;
    for K = 1 : length(Aval); hidden_h(K) = surf(uint8(K-[1 1;1 1]), 'edgecolor', 'none'); end
    hold off
    uistack(hidden_h, 'bottom');
    legend(hidden_h, tissues,'NumColumns',2)
    clear hidden_h
    title('Original segmented image')
    xlabel('cm')
    ylabel('cm')
    axis on
    xticks(linspace(1,512,length(0:10:51.2)))
    xticklabels(0:10:51.2)
    yticks(linspace(1,512,length(0:10:51.2)))
    yticklabels(0:10:51.2)

    % create the phantom image by replacing the ids with mean values of regions
    % of the transformed attenuation coefficient matrix 
    P = zeros(size(P_ct));
    for i = 1:length(Aval)
        P(Anew==i) = mean(P_ct_mu(Anew==i));
    end
    subplot(1,3,3)
    imshow(P,[])
    title('Transformed phantom image')
    xlabel('cm')
    ylabel('cm')
    axis on
    xticks(linspace(1,512,length(0:10:51.2)))
    xticklabels(0:10:51.2)
    yticks(linspace(1,512,length(0:10:51.2)))
    yticklabels(0:10:51.2)
    subplot(1,3,1)
    imshow(P_ct_mu,[])
    title('Original reconstructed image')
    xlabel('cm')
    ylabel('cm')
    axis on
    xticks(linspace(1,512,length(0:10:51.2)))
    xticklabels(0:10:51.2)
    yticks(linspace(1,512,length(0:10:51.2)))
    yticklabels(0:10:51.2)
    
    % scale the phantom for (mm-1) unit
    P = P/10;
end

% create the attenuation coefficient transform polynomials from abdomen CT
function [p] = mu_polyfit
    % choose a slice which has air, kidney and muscle regions
    segment = 'BODY';
    slice = 40;

    % read CT image file convert
    fileID = fopen(strcat('dataset\CT_slices\',segment,'\',segment,'_',num2str(slice,'%03d'),'_I.DAT.1'));
    P_ct = fread(fileID,[512,512],'int16=>int16');
    P_ct = double(swapbytes(P_ct)).';

    % read color segmented image file convert
    fileID = fopen(strcat('dataset\Color_slices\',segment,'\',segment,'_',num2str(slice,'%03d'),'_C.DAT.1'));
    P_clr = fread(fileID,[512,512],'uint8=>double').';

    % possible efficient energy levels
    keVs = [60,65,70,75,80,85,90];
    mu_air = [0,0,0,0,0,0,0];
    mu_kidney = [0.2127,0.2071,0.2013,0.1958,0.1908,0.1875,0.1842]; % these values are from a paper
    mu_muscle = [0.2111,0.2052,0.1993,0.1935,0.1878,0.1844,0.1812]; % these values are from a paper
    p = zeros(length(keVs),2); % fit polynomials
    s = zeros(1,length(keVs)); % fit errors
    for i = 1:length(keVs) 
        % find the mean values in CT image of the regions of segmented image 
        x = [mean(P_ct(P_clr==0)),mean(P_ct(P_clr==77)),mean(P_ct(P_clr==72))];
        y = [mu_air(i),mu_kidney(i),mu_muscle(i)];
        [p(i,:),S] = polyfit(x,y,1); % find the best fit
        s(i) = S.normr/mean(y);
        fprintf("%d keV, error: %d\n",keVs(i),s(i))
    end 
    % the best fit with the least error is the efficient energy
    [~,I] = min(s);
    fprintf("\nThe best fit is at %d keV\n",keVs(I))
    % the corresponding polynomials are used in the transform 
    p = p(I,:);
end