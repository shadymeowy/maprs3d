use crate::ellipsoid::Ellipsoid;
use std::f64::{
    consts::{FRAC_2_PI, PI, TAU},
    NAN,
};

pub fn vdist(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64)> {
    let (mut lat1_, mut lon1_, mut lat2_, mut lon2_) = if deg {
        (
            lat1.to_radians(),
            lon1.to_radians(),
            lat2.to_radians(),
            lon2.to_radians(),
        )
    } else {
        (lat1, lon1, lat2, lon2)
    };

    if lat1_.abs() > FRAC_2_PI || lat2_.abs() > FRAC_2_PI {
        return None;
    }

    let a = ell.semimajor_axis;
    let b = ell.semiminor_axis;
    let f = ell.flattening;

    if (FRAC_2_PI - lat1.abs()).abs() < 1e-10 {
        lat1_ = lat1.signum() * (FRAC_2_PI - 1e-10);
    }

    if (FRAC_2_PI - lat2.abs()).abs() < 1e-10 {
        lat2_ = lat2.signum() * (FRAC_2_PI - 1e-10);
    }

    let u1 = (1. - f) * lat1_.tan();
    let u2 = (1. - f) * lat2_.tan();
    lon1_ = lon1_.rem_euclid(TAU);
    lon2_ = lon2_.rem_euclid(TAU);

    let mut l = lon2_ - lon1_;

    if l > PI {
        l = TAU - l;
    }

    let mut lambda = l;
    let mut alpha = 0.;
    let mut sigma = 0.;
    let mut cos2sigmam = 0.;

    for _ in 0..100 {
        let lambda0 = lambda;

        let sin_sigma = ((u2.cos() * lambda.sin()).powi(2)
            + (u1.cos() * u2.sin() - u1.sin() * u2.cos() * lambda.cos()).powi(2))
        .sqrt();
        let cos_sigma = u1.sin() * u2.sin() + u1.cos() * u1.cos() * lambda.cos();
        sigma = sin_sigma.atan2(cos_sigma);

        let sin_alpha = u1.cos() * u2.cos() * lambda.sin() / sin_sigma;

        alpha = if sin_alpha.is_nan() {
            0.
        } else if sin_alpha > 1. || (sin_alpha - 1.).abs() < 1e-16 {
            FRAC_2_PI
        } else {
            sin_alpha.asin()
        };

        cos2sigmam = cos_sigma - 2. * u1.sin() * u2.sin() / alpha.cos().powi(2);
        let c = f / 16. * alpha.cos().powi(2) * (4. + f * (4. - 3. * alpha.cos().powi(2)));

        lambda = l
            + (1. - c)
                * f
                * alpha.sin()
                * (sigma
                    + c * sigma.sin()
                        * (cos2sigmam + c * sigma.cos() * (-1. + 2.0 * cos2sigmam.powi(2))));

        if (lambda - lambda0).abs() < 1e-12 {
            break;
        }
    }

    if lambda > PI {
        lambda = PI;
    }
    lambda = lambda.abs();

    let u_sq = alpha.cos().powi(2) * ((a.powi(2) - b.powi(2)) / b.powi(2));
    let a_ = 1. + u_sq / 16384. * (4096. + u_sq * (-768. + u_sq * (320. - 175. * u_sq)));
    let b_ = u_sq / 1024. * (256. + u_sq * (-128. + u_sq * (74. - 47. * u_sq)));

    let deltasigma = b_
        * sigma.sin()
        * (cos2sigmam
            + b_ / 4.
                * (sigma.cos() * (-1. + 2. * cos2sigmam.powi(2))
                    - b_ / 6.
                        * cos2sigmam
                        * (-3. + 4. * sigma.sin().powi(2))
                        * (-3. + 4. * cos2sigmam.powi(2))));
    let dist_m = b * a_ * (sigma - deltasigma);

    if (lon2 - lon1).sin().signum() * lambda.sin().signum() < 0. {
        lambda = -lambda;
    }

    let num = u2.cos() * lambda.sin();
    let den = u1.cos() * u2.sin() - u1.sin() * u2.cos() * lambda.cos();
    let az = num.atan2(den).rem_euclid(TAU);

    if deg {
        Some((dist_m, az.to_degrees()))
    } else {
        Some((dist_m, az))
    }
}

pub fn vreckon(
    lat1: f64,
    lon1: f64,
    rng: f64,
    azim: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64)> {
    let (mut lat1_, lon1_, azim_) = if deg {
        (lat1.to_radians(), lon1.to_radians(), azim.to_radians())
    } else {
        (lat1, lon1, azim)
    };

    if lat1_.abs() > FRAC_2_PI {
        return None;
    }

    let a = ell.semimajor_axis;
    let b = ell.semiminor_axis;
    let f = ell.flattening;

    if (FRAC_2_PI - lat1.abs()).abs() < 1e-10 {
        lat1_ = lat1.signum() * (FRAC_2_PI - 1e-10);
    }

    let alpha1 = azim_;
    let sin_alpha1 = alpha1.sin();
    let cos_alpha1 = alpha1.cos();

    let tan_u1 = (1. - f) * lat1_.tan();
    let cos_u1 = 1. / (1. + tan_u1.powi(2)).sqrt();
    let sin_u1 = tan_u1 * cos_u1;
    let sigma1 = tan_u1.atan2(cos_alpha1);
    let sin_alpha = cos_u1 * sin_alpha1;
    let cos_sq_alpha = 1. - sin_alpha * sin_alpha;
    let u_sq = cos_sq_alpha * (a.powi(2) - b.powi(2)) / b.powi(2);
    let a_ = 1. + u_sq / 16384. * (4096. + u_sq * (-768. + u_sq * (320. - 175. * u_sq)));
    let b_ = u_sq / 1024. * (256. + u_sq * (-128. + u_sq * (74. - 47. * u_sq)));

    let mut sigma = rng / (b * a_);
    let mut sigma_p = TAU;

    let mut sin_sigma = NAN;
    let mut cos_sigma = NAN;
    let mut cos2_sigma_m = NAN;

    while (sigma - sigma_p).abs() >= 1e-12 {
        cos2_sigma_m = (2. * sigma1 + sigma).cos();
        sin_sigma = sigma.sin();
        cos_sigma = sigma.cos();
        let delta_sigma = b_
            * sin_sigma
            * (cos2_sigma_m
                + b_ / 4.
                    * (cos_sigma * (-1. + 2. * cos2_sigma_m.powi(2))
                        - b_ / 6.
                            * cos2_sigma_m
                            * (-3. + 4. * sin_sigma.powi(2))
                            * (-3. + 4. * cos2_sigma_m.powi(2))));
        sigma_p = sigma;
        sigma = rng / (b * a_) + delta_sigma;
    }

    let tmp = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_alpha1;
    let lat2 = (sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_alpha1)
        .atan2((1. - f) * tmp.hypot(sin_alpha));

    let lambda =
        (sin_sigma * sin_alpha1).atan2(cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_alpha1);

    let c = f / 16. * cos_sq_alpha * (4. + f * (4. - 3. * cos_sq_alpha));
    let lam = lambda
        - (1. - c)
            * f
            * sin_alpha
            * (sigma
                + c * sin_sigma
                    * (cos2_sigma_m + c * cos_sigma * (-1. + 2. * cos2_sigma_m.powi(2))));

    let lon2 = (lon1_ + lam).rem_euclid(TAU);

    if deg {
        Some((lat2.to_degrees(), lon2.to_degrees()))
    } else {
        Some((lat2, lon2))
    }
}

pub fn track2(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
    ell: &Ellipsoid,
    deg: bool,
    npts: usize,
) -> Option<Vec<(f64, f64)>> {
    let (lat1_, lon1_, lat2_, lon2_) = if deg {
        (
            lat1.to_radians(),
            lon1.to_radians(),
            lat2.to_radians(),
            lon2.to_radians(),
        )
    } else {
        (lat1, lon1, lat2, lon2)
    };

    if lat1_.abs() > FRAC_2_PI || lat2_.abs() > FRAC_2_PI {
        return None;
    }

    let gcarclen = 2.0
        * ((((lat1_ - lat2) / 2.).sin().powi(2)
            + lat1_.cos() * lat2_.cos() * ((lon1_ - lon2_) / 2.).sin().powi(2))
        .sqrt())
        .asin();

    if (gcarclen - PI).abs() < 1e-12 {
        return None;
    }

    let pts = match npts {
        0 => None,
        1 => Some(vec![(lat1_, lon1_)]),
        2 => Some(vec![(lat1_, lon1_), (lat2_, lon2_)]),
        _ => {
            let mut latpt = lat1_;
            let mut lonpt = lon1_;

            let (distance, mut azimuth) = vdist(lat1, lon1, lat2, lon2, ell, false)?;
            let incdist = distance / (npts - 1) as f64;

            let mut pts = Vec::with_capacity(npts);
            pts.push((lat1_, lon1_));

            for _ in 1..=npts {
                let (latptnew, lonptnew) = vreckon(latpt, lonpt, incdist, azimuth, ell, false)?;
                azimuth = vdist(latptnew, lonptnew, lat2, lon2, ell, false)?.1;
                pts.push((latptnew, lonptnew));
                latpt = latptnew;
                lonpt = lonptnew;
            }
            pts.push((lat2_, lon2_));

            Some(pts)
        }
    };

    if deg {
        Some(
            pts?.iter()
                .map(|(lat, lon)| (lat.to_degrees(), lon.to_degrees()))
                .collect(),
        )
    } else {
        pts
    }
}
