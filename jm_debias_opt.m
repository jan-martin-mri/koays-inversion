function opt = jm_debias_opt( opt )
% jm_debias_opt Generate options for jm_debias.
%     do_noisemap: Generate additional noise map during debiasing
%     do_snrmap: Generate additional snr map during debiasing
%     fn_debias: File name attachment for debiased data
%     fn_noisemap: File name attachment for generated noise map
%     fn_snrmap: File name attachment for generated snr map
%     ref_values: File name for lookup table containing the corrected snr
%     values
opt.present = 1;

%% Debiasing
opt.debias.present = 1;
opt.debias = msf_ensure_field(opt.debias, 'do_noisemap', 0); % 0, 1
opt.debias = msf_ensure_field(opt.debias, 'do_snrmap', 0); % 0, 1
opt.debias = msf_ensure_field(opt.debias, 'fn_debias', '_db');
opt.debias = msf_ensure_field(opt.debias, 'fn_noisemap', '_db_noisemap');
opt.debias = msf_ensure_field(opt.debias, 'fn_snrmap', '_db_snr');
opt.debias = msf_ensure_field(opt.debias, 'ref_values', 'jm_debias_ref.mat');

end

