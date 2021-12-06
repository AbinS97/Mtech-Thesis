function [first_ht] =  multi_ht_iteration(N, wavelength, area_size,z,heights, iterations )

%N : pixel numbers along a dimension 
%wavelength: wavelength of light used(in m)
%area_size : the side length of detector(in metres)
%z : the distance between the hologram recordings
%heights : the number of heights taken 
% iterations : the number of iterations doing along the holograms
% Z : the distance between fist height and the object plane. (Not needed here)

for hh = 1: heights
    
    fid = fopen(['a',num2str(hh),'_hologram.bin'], 'r');
    hologram = fread(fid, [N, N], 'real*4');
    fclose(fid); 
    measured{hh}= sqrt(hologram);
end

phase = zeros(N,N); %initial assumption of phase value
first_ht = measured{1}.*exp(i*phase); %initial assumption at height 1
prop = Propagator(N, wavelength, area_size, z);

for tt=1: iterations
    disp(tt)
for nn = 1: (heights -1)
    
t = IFT2Dc(FT2Dc(first_ht).*conj(prop));
%imshow(angle(t))
angle(t);
am = (measured{nn+1}+abs(t))/2;
ph = angle(t);
first_ht = am.* exp(i*ph);
end
%imshow
for mm = 1: (heights -1)
t = IFT2Dc(FT2Dc(first_ht).*(prop));
am =( measured{heights - mm}+abs(t))/2;
ph = angle(t);
%imshow(ph)
%imshow(ph)
first_ht = am.* exp(i*ph);
 
end
end
end

