%% simulating image formation using a microscope

% Alex Kariuki Muthumbi, 3rd September, 2019

%% Steps to take
%   - create the spatial coordinates, illuminating beam and simulate image being
%     illuminated by the beam. 
%   - load in the images to use for simulating.
%   - Create Fourier coordinates, and take the FFT of the illuminated
%     sample to take the field to the lens
%   - Create the transfer function of the lens and simulate the effect of
%   the lens on the field at the lens
%   - Take Inverse FFT to take the field to the image plane
%   - Simulate the sampling of the field by the camera


%% load in the images to use for simulating.
% i will 2 images readily available in matlab

amp = double(imread('cameraman.tif'));

% use pepper as phase image
phase = double(imread('peppers.png')); % this is a RGB image.
% Use the red channel only
phase = phase(:,:,1);

% let both images to have same dimensions 
phase = imresize(phase, [size(amp,1) size(amp,2)]);

%% Define spatial coordinates
wavelength = 0.5e-3 ; % wavelength to be used, in mm
k0 = 2*pi / wavelength; % define wave vector as its used in many places

%  sample step
sample_step = wavelength / 2; % defined from the sampling theory

% define number of samples to be used. Should be same size as the sample
num_of_samples = size(amp, 1);

% define spatial coordinates
starting_coordinate = (-num_of_samples/2) * sample_step; 
last_coordinate = (num_of_samples/2) * sample_step;

% create equally spaced samples
x = linspace(starting_coordinate, last_coordinate, num_of_samples); 
y = linspace(starting_coordinate, last_coordinate, num_of_samples);

% define a meshgrid of the samples
[xx, yy] = meshgrid(x,y);

% define the angle in degrees of the illuminating wave, and convert it to radians
angle_x_in_degrees = 0; % can change this if you want illumination from certain angle
angle_x_in_radians = angle_x_in_degrees / pi * 180;
angle_y_in_degrees = 0; % can change this if you want illumination from certain angle
angle_y_in_radians = angle_y_in_degrees / pi * 180; 
% create a plane wave to illuminate the sample
illuminating_wave = exp(1i * k0 * ((sin(angle_x_in_radians) * xx) + sin(angle_y_in_radians) * yy));

% define constant alpha to be used in making sure that there is no phase
% wrapping
alpha = 1e-8;

% create a complex image now
sample = amp .* exp(1i * k0 * alpha * phase); 

% see the complex sample now
figure; imagesc(abs(sample)); axis image; colormap(gray(256));title('Amplitude of sample')
figure; imagesc(angle(sample)); axis image; colormap(gray(256));title('Phase of sample')

% Simulate illumination of the sample with the illuminating wave
illuminated_sample = sample .* illuminating_wave;

% figure; imagesc(abs(illuminated_sample)); axis image; colormap(gray(256)); title('Amplitude of illuminated sample')
% figure; imagesc(angle(illuminated_sample)); axis image; colormap(gray(256)); title('Phase of illuminated sample')

% see the illuminating angle
figure; imagesc(angle(illuminating_wave)); axis image; colormap(gray(256));

%% Fourier plane

% define the Fourier coordinates

% The range in Fourier plane is equal to the inverse of the step in spatial
% coordinates (sample_step)
fx_range = 1/sample_step ;  % dimensions in mm^-1

% the step in Fourier plane is inverse of the total range in spatial
% coordinates (sample_step * num_of_samples)
dfx = 1/ (sample_step * num_of_samples); % dimensions also in mm^-1

% define now the fx
starting_fx = -fx_range/2;
last_fx = fx_range/2;

fx = linspace(starting_fx, last_fx, num_of_samples);  %  fx = fy for us

% Take the Fourier transform of the sample to move from the illuminating
% plane to the Fourier plane
F_sample = fftshift(fft2(fftshift(illuminated_sample)));
figure; imagesc(fx, fx, log(abs(F_sample))); axis image; colormap(gray(256)); title('Amp of Fourier');xlabel('mm ^{-1}'); ylabel('mm ^{-1}');

% we can see the Fourier plane in 2 ways, 1 in fx coordinates, or 2 in
% angular coordinates. 
% angular coordinates are connected to the Fourier coordinates using the
% equation fx = sin(angle_theta)/ wavelength
% therefore, we can get the sin(angle_theta) by

sin_theta = fx .* wavelength;

% figure; imagesc(sin_theta, sin_theta, log(abs(F_sample))); axis image; colormap(gray(256)); title('Amp of Fourier');xlabel('sin \theta'); ylabel('sin \theta');


%% Lens effect

% The lens low pass filters the field at the Fourier plane. To define how
% this works, we need to define the transfer function of the lens, which is
% dependent on the lens' Numerical apereture 

lens_NA = 0.4;

% define the transfer function
[sin_theta_x, sin_theta_y] = meshgrid(sin_theta, sin_theta);
lens_transfer_function = double(sqrt(sin_theta_x.^2 + sin_theta_y.^2) < lens_NA);
figure; imagesc(lens_transfer_function); axis image; colormap(gray(256)); title('Lens transfer function')


% low pass filter the sample now
F_sample_filtered = F_sample .* lens_transfer_function;
figure; imagesc(log(abs(F_sample_filtered))); axis image; colormap(gray(256)); title('Filtered Fourier transform of image')

% We may also want to simulate the effect of magnification of the field by
% the lens as some lenses have a magnification greater than 1

lens_magnification = 2;

% magnify the field now
F_sample_magnified = imresize(F_sample_filtered, [size(F_sample_filtered,1)*lens_magnification size(F_sample_filtered,2)*lens_magnification]);
figure; imagesc(log(abs(F_sample_magnified))); axis image; colormap(gray(256)); title('Magnified filtered Fourier transform of image')

%% Image plane

% Take the inverse Fourier transform of the filtered Fourier transform of
% the sample to move to the 
image_plane_field = ifftshift(ifft2(ifftshift(F_sample_magnified))); 

% At this plane, we normally have a CCD or CMOS camera. These types of
% cameras can only detect the field's intenstity field, which is the
% absolute value of the field at the plane squared
image_intensity = abs(image_plane_field).^2;

figure; imagesc(image_intensity); axis image; colormap(gray(256)); title('Intensity at ')

%% Sampling by camera

% the camera discretely samples the field at the image plane, depending on
% the number of pixels on the camera

pixel_x = 1024;
pixel_y = 768;

% preallocate, to save on memory
final_camera_image = zeros(pixel_x, pixel_y);

% create equally spaced values 
[detector_x] = linspace(1, pixel_x, pixel_x); % 1 to 1024, 1024 points
[detector_y] = linspace(1, pixel_y, pixel_y);

% resize the dimensions above to the same size as the intensity image
resized_pixel_x = imresize(detector_x, [1 size(image_intensity, 1)], 'nearest');
resized_pixel_y = imresize(detector_y, [1 size(image_intensity, 1)], 'nearest');

% resample the field now
for aa = 1:size(image_intensity, 1)
    for bb = 1: size(image_intensity, 2)
        
        final_camera_image(resized_pixel_x(aa), resized_pixel_y(bb)) = final_camera_image(resized_pixel_x(aa), resized_pixel_y(bb)) + image_intensity(aa,bb);
        
    end 
end

% view the image recorded by camera
figure; imagesc(final_camera_image); axis image; colormap(gray(256)); title('Image recorded by camera')
