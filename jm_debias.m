function jm_debias(fn_data, fn_noisemap, opt)
% jm_debias Correct data for Rician noise distribution using Koay's
% inversion method.
%     This function applies Koay's inversion method to noisy data using a
%     pre-calculated look-up table from an external file to speed up the 
%     process (see end of file for calculation details).
% Input:
%     fn_data: Path to nii.gz file containing noisy data
%     fn_noisemap: Path to nii.gz file containing noisemap corresponding to
%     noisy data
%     opt (optional input): Options struct containing processing options
%     for jm_debias (see jm_debias_opt for more detail)
% Output files are generated in the same folder as the input noisy data.

% Processing options
if nargin < 3
    opt = jm_debias_opt();
else
    opt = jm_debias_opt( opt );
end

% Load precalculated reference SNR values
ref_input_snr = [];
ref_output_snr = [];
load(opt.debias.ref_values, 'ref_input_snr', 'ref_output_snr');

% Load signal and noise data
[input_signal, h] = mdm_nii_read(fn_data);
[input_noise, h_noise] = mdm_nii_read(fn_noisemap);
    % Rescale noise data to match signal data
if size(input_noise,4) == 1
    input_noise = repmat(input_noise,1,1,1,size(input_signal,4));
end

% Check that the dimensions of the signal data, noise data, and mask are
% matching
if ~isequal(size(input_signal), size(input_noise))
    error('Input signal and noise map do not have matching dimensions.');
end

% Save the dimensions of the signal data to rescale the results later
dim = size(input_signal);

input_signal = double(input_signal(:));
input_noise  = double(input_noise(:));
input_snr = input_signal ./ input_noise; 

output_signal = zeros(size(input_signal));
output_noise  = ones(size(input_noise));
output_snr    = zeros(size(input_snr));

% Separate the SNR into three distinct regions:
%   a. SNR < 1.913 which can not be corrected and will subsequently be set
%      to 0
%   b. SNR > 1.913 but < 10 which will be corrected using the look-up
%      tables in jm_debias_ref
%   c. SNR > 10 or SNR = NaN which will not be corrected
ind_a = input_snr < min(ref_input_snr);
ind_b = (input_snr >= min(ref_input_snr)) & (input_snr <= max(ref_input_snr));
ind_c = (input_snr > max(ref_input_snr)) | isnan(input_snr);

% Check if all elements of input_snr are accounted for
if sum(ind_a + ind_b + ind_c) ~= numel(input_snr)
    error('Indexing of input_snr went wrong. Please check jm_debias for possible bugs');
end

% b. Correct SNR using the look-up tables in jm_debias_ref
snr_b = input_snr(ind_b);
% Ensure that snr_b and ref_input_snr can be compared
snr_b = round(snr_b, 3);
ref_input_snr = round(ref_input_snr, 3);
% Reduce snr_b to its unique entries
[snr_b, ~, ind_b_unique] = unique(snr_b);
% Replace the original SNR values with the corrected SNR values
snr_b = ref_output_snr(ismember(ref_input_snr, snr_b));
% Calculate the signal based on the corrected SNR values
e = @(x) 2 + x.^2 - pi / 2 * (hypergeom(-1/2, 1, -(x.^2)/2)).^2;
e_snr_b = e(snr_b');
% Rescale e_snr_b
% Careful: Unique SNR values can still be composed of several different
% combinations of signal and noise values. Therefore signal and noise may
% not be reduced to unique values in the same way as the SNR.
e_snr_b = e_snr_b( ind_b_unique );
% Extract signal and noise values for interval b.
input_signal_b = input_signal( ind_b );
input_noise_b  = input_noise( ind_b );
% Calculate the corrected signal and noise
signal_b = sqrt(input_signal_b.^2 + (1 - 2 ./ e_snr_b) .* input_noise_b.^2);
if opt.debias.do_noisemap
    noise_b = sqrt(input_noise_b.^2 ./ e_snr_b);
end

% Export the corrected signal
output_signal(ind_a) = 0;
output_signal(ind_b) = signal_b;
output_signal(ind_c) = input_signal(ind_c);
output_signal = reshape(output_signal, dim);

[pn, fn, ext] = msf_fileparts(fn_data);
fn_debias = fullfile(pn, strcat(fn, opt.debias.fn_debias, ext));
mdm_nii_write(output_signal, fn_debias, h);

% Export the corrected SNR
if opt.debias.do_snrmap
    snr_b = snr_b( ind_b_unique );
    output_snr(ind_a) = 0;
    output_snr(ind_b) = snr_b;
    output_snr(ind_c) = input_snr(ind_c);
    output_snr = reshape(output_snr, dim);
    
    fn_snrmap = fullfile(pn, strcat(fn, opt.debias.fn_snrmap, ext));
    mdm_nii_write(output_snr, fn_snrmap, h);
end

% Export the corrected noise
if opt.debias.do_noisemap
    noise_b = noise_b( ind_b_unique );
    output_noise(ind_a) = 1;
    output_noise(ind_b) = noise_b;
    output_noise(ind_c) = input_noise(ind_c);
    output_noise  = reshape(output_noise, dim);

    fn_noisemap = fullfile(pn, strcat(fn, opt.debias.fn_noisemap, ext));
    mdm_nii_write(output_noise, fn_noisemap, h_noise);
end

end

%% Debiased SNR values were calculated according to Koay's inversion:
% input_snr = 1.914 : 0.001 : 100;
% 
% % Define functions for root-finding algorithm (see Koay, 2006)
% e = @(x) 2 + x^2 - pi / 2 * (hypergeom(-1/2, 1, -(x^2)/2))^2;
% g = @(x,r) sqrt(e(x) * (1 + r^2) - 2);
% 
% output_snr = zeros(size(input_snr));
% 
% ppm = ParforProgressbar(numel(input_snr), 'progressBarUpdatePeriod', 1);
% options = optimoptions('fsolve', 'UseParallel', true, 'Display', 'off');
% 
% parfor idx = 1 : numel(input_snr)
%     % Only run the algorithm if the SNR is higher than the
%     % lowerBound of 1.9130
%     if (input_snr(idx) > 1.913)
%         r = input_snr(idx);
%         fun = @(x) g(x,r) - x;
%         output_snr(idx) = fsolve(fun, r - 1.9130, options);
%     end
%     % Update parfor progress monitor
%     ppm.increment();
% end
% 
% % Clean-Up parfor progress monitor
% delete(ppm);