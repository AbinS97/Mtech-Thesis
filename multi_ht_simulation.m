
% Get 5 different heights of the hologram 

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters

N = 500;                    % number of pixels
lambda = 532*10^(-9);       % wavelength in meter
object_area = 0.002;        % object area sidelength in meter
z = 0.005; %1mm         %0.05         % object-to-detector distance in meter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating object pxccvvlane transmission function

object = zeros(N,N);        
    object0 = imread('a_hair.jpg');   
    object(:,:) = object0(:,:,1);    
    object = (object - min(min(object)))/(max(max(object)) - min(min(object)));    
    figure, imshow(object);  
    
    


am = exp(-1.6*object);    % exp(-1.6)=0.2 %why is it multiplied by -1.6? 
figure, imshow(am)
ph = - 3*object; %is it assumed phase? Yes 

t = zeros(N,N);

t = am.*exp(i*ph); 
figure, imshow(angle(t))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating hologram
 

iterations = 8; 
for kk = 1: iterations 
    prop = Propagator(N, lambda, object_area, (z + ((kk-1)*90*(10^-6))));
    U = IFT2Dc(FT2Dc(t).*conj(prop));
    hologram = abs(U).^2;

    figure('Name',['Simulated hologram','_height_',num2str(kk)],'NumberTitle','off')
    imshow(hologram, []), colorbar;
    colormap(gray)  
    fid = fopen(strcat(['a',num2str(kk),'_hologram.bin']), 'w');
    fwrite(fid, hologram, 'real*4');
    fclose(fid);


    p = hologram;
    p = 255*(p - min(min(p)))/(max(max(p)) - min(min(p)));
    imwrite (p, gray, ['a',num2str(kk),'_hologram.bmp']);

    
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving hologram


