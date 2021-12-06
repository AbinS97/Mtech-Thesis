

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 500;                    % number of pixels
Iterations = 30;            % number of iterations

wavelength = 532*10^(-9);   % wavelength in meter
area_size = 0.002;          % area size = object area size, in meter
                            % object-to-detector distance in meter

p = 0.01;                   % time to pause, otherwise images are not shown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reading all hologram intensities to the variables at the 3 heights
heights = 8;
for hh = 1: heights
    
    fid = fopen(['a',num2str(hh),'_hologram.bin'], 'r');
    hologram = fread(fid, [N, N], 'real*4');
    fclose(fid); 
    measured{hh}= sqrt(hologram)
end

% Creating initial complex-valued field distribution in the detector plane
phase = zeros(N,N); %initial assumption of phase value

frst_ht = measured{1}.*exp(i*phase); %initial assumption at height 1

% Creating wave propagation term
%prop = Propagator(N, wavelength, area_size, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative reconstruction
 
%figure('Name','Current amplitude = exp(-absorption) / a.u.','NumberTitle','off')
 
%iterate through all the heights 
z = 30*(10^-6);
prop = zeros(N,N);

ht_iterations = 30;

for kk = 1:ht_iterations

%the count of iterations
fprintf('Iteration: %d\n', kk) 
%field_detector = measured.*exp(i*phase);

% forward propagating to the 2nd height

prop = Propagator(N, wavelength, area_size, z);
t = IFT2Dc(FT2Dc(frst_ht).*conj(prop));


%updating the value at the second height 

am = measured{2};
ph = angle(t);
scnd_ht  = am.*exp(i*ph);

%finding the value at the 3rd height 
t = IFT2Dc(FT2Dc(scnd_ht).*conj(prop));


%updating the amplitude

am = measured{3};
ph = angle(t);
thrd_ht = am.*exp(i*ph);



%at 4th height 

t = IFT2Dc(FT2Dc(thrd_ht).*conj(prop));

am = measured{4};
ph = angle(t);
frth_ht = am.*exp(i*ph);

%5th height 
t = IFT2Dc(FT2Dc(frth_ht).*conj(prop));

am = measured{5};
ph = angle(t);
fifth_ht = am.*exp(i*ph);

%back propagation
%back propagating to 4th height 


t = IFT2Dc(FT2Dc(fifth_ht).*prop);
am = measured{4};
ph = angle(t);
frth_ht = am.*exp(i*ph);

%back propagating to 3rd height 


t = IFT2Dc(FT2Dc(frth_ht).*prop);
am = measured{3};
ph = angle(t);
thrd_ht = am.*exp(i*ph);

%back propagating to 2nd height 
t = IFT2Dc(FT2Dc(thrd_ht).*prop);
am = measured{2};
ph = angle(t);
scnd_ht = am.*exp(i*ph);

%back propagating to 1st height

t = IFT2Dc(FT2Dc(scnd_ht).*prop);
am = measured{1};
ph = angle(t);
frst_ht = am.*exp(i*ph);



end

%updating the information at the 2nd height 






%imshow(am, [],'colormap', 1-gray)
%pause(p);

%backpropagating to the object plane 
z = 0.05  ; 
prop = Propagator(N, wavelength, area_size, z);
%first height to object distance 
object_plane = IFT2Dc((FT2Dc(frst_ht)).*prop);

amplitude = abs(object_plane);
phase = angle(object_plane);

figure('Name','first height amplitude ')
imshow(am);

figure('Name','first height phase ')
imshow(ph);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Showing reconsrtucted amplitude
    figure('Name','Reconstructed amplitude = exp(-absorption) / a.u.','NumberTitle','off')
    imshow(amplitude, [],'colormap', 1-gray), colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Showing reconsrtucted phase
    figure('Name','Reconstructed phase / radian','NumberTitle','off')
    imshow(phase, [],'colormap', 1-gray), colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving reconsrtucted amplitude as jpg image
        h = - amplitude;
        h = (h - min(min(h)))/(max(max(h)) - min(min(h)));
        imwrite (h, 'reconstructed_amplitude_multi_ht.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving reconsrtucted phase as jpg image
        h = - phase;
        h = (h - min(min(h)))/(max(max(h)) - min(min(h)));
        imwrite (h, 'reconstructed_phase_multi_ht.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
