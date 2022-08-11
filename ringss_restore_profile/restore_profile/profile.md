### Cube parameters

```
cube = { 
            'nx': roi_wdth_int,
            'nz': number_of_cubes_int,
            'texp': exposure_usec_int,
            'gain': cmos_acq_gain_int,
            'data': '2022-03_26T19:30:00',
            'star': HR_int
            }
```
### Star parameters
```
star = {
            'JD': julian_time_float,
            'HR': HR_int,
            'Vmag': visual_magnitude_float,
            'BV': B-V_color_float,
            'zen': zenith_deg_float,
            'az': azimuth_deg_float
            }
```         

### Image parameter

```
# All values are float

Image_parameter = {
            'image': 
                    {
                    'impar': [backgr, flux, fluxvar, meanrad, rwidth, xc, yc, xcvar, ycvar, coma, angle, contrast], 
                    'noisepar': [noise1, noise2, noise1r, noise2r]
                    }, 
            'moments': 
                    {
                    'var': [list_of_variances],
                    'cov': [list_of_covariances]
                    'rnoise': signal_noise_float,
                    'rvar': mean_ring_var_float,
                    'mcoef': [list_of_mean_coef_for_aberrations] 
                    }
            }


```
### Restore Profile input data
An extended version of Image_parameter dictionary that includes cube parameters and star parameters. 
```
Data = {
        'image' : { },
        'moments': { },
        'cubepar' { },
        'star_par': { }
}

```

### Turbulence Profile output
```
# All values are float

Turbulence_profile = 
            {
            "z0": z0, 
            "prof": [list], 
            "see": total_seeing_arcsecs, 
            "fsee": fa_seeing_arcsecs, 
            "see2": see2, 
            "wind": v2,
            "erms": erms,
            "totvar": totvar, 
            "tau0": tau0, 
            "theta0": theta0
            }

```
## Profile Restoration

![my profile restoration schema](./RINGSS_profile_restoration_offline_diagram.png "RINGSS data path")
![my profile restoration schema](./RINGSS_profile_restoration_online_diagram.png "RINGSS data path")
