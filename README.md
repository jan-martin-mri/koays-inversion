# Koay's Inversion Method

This GitHub repository contains MATLAB code for correcting noisy MRI data for Rician bias according to Koay's inversion method (https://doi.org/10.1016/j.jmr.2006.01.016).

Using Koay's inversion method to correct data for Rician bias requires a corresponding noise map. Such a noise map may for example be generated using freely available denoising methods such as "dwidenoise" which is part of the MRtrix3 v3.0.2 toolbox (https://www.mrtrix.org/download/). The data and the noise map need to be formatted as .nii.gz.

Example usage:

```matlab
opt = jm_debias_opt();
jm_debias(fn_data, fn_noisemap, opt);
```

Using this toolbox requires the "md-dmri" toolbox available at https://github.com/jan-martin-mri/md-dmri.