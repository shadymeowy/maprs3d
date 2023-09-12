use crate::ellipsoid::Ellipsoid;

pub fn geodetic2spherical(
    lat: f64,
    lon: f64,
    alt: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let lat = if deg { lat.to_radians() } else { lat };
    let lon = if deg { lon.to_radians() } else { lon };

    if lat.abs() > std::f64::consts::FRAC_PI_2 {
        return None;
    }

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();

    let n = ell.semimajor_axis.powi(2)
        / (ell.semimajor_axis * cos_lat).hypot(ell.semiminor_axis * sin_lat);

    let xy = (n + alt) * cos_lat;
    let z = (alt + (1. - ell.eccentricity.powi(2)) * n) * sin_lat;
    let r = xy.hypot(z);
    let slat = (z / r).asin();

    let slat = if deg { slat.to_degrees() } else { slat };
    let lon = if deg { lon.to_degrees() } else { lon };

    Some((slat, lon, r))
}

pub fn spherical2geodetic(
    lat: f64,
    lon: f64,
    r: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let lat = if deg { lat.to_radians() } else { lat };
    let lon = if deg { lon.to_radians() } else { lon };

    if lat.abs() > std::f64::consts::FRAC_PI_2 {
        return None;
    }

    let sin_lat = lat.sin();
    let cos_lat = lat.cos();

    let z = r * sin_lat;
    let p_0 = r.powi(2) * cos_lat.powi(2) / ell.semimajor_axis.powi(2);
    let q_0 = (1. - ell.eccentricity.powi(2)) / ell.semimajor_axis.powi(2) * z.powi(2);
    let r_0 = (p_0 + q_0 - ell.eccentricity.powi(4)) / 6.;
    let s_0 = ell.eccentricity.powi(4) * p_0 * q_0 / 4. / r_0.powi(3);
    let t_0 = 1. + s_0 + (2. * s_0 + s_0.powi(2)).cbrt().sqrt();
    let u_0 = r_0 * (1. + t_0 + 1. / t_0);
    let v_0 = (u_0.powi(2) + q_0 * ell.eccentricity.powi(4)).sqrt();
    let w_0 = ell.eccentricity.powi(2) * (u_0 + v_0 - q_0) / 2. / v_0;
    let k = (u_0 + v_0 + w_0.powi(2)).sqrt() - w_0;
    let d = k * r * cos_lat / (k + ell.eccentricity.powi(2));
    let hypot_dz = d.hypot(z);

    let glat = 2. * z.atan2(d + hypot_dz);
    let alt = (k + ell.eccentricity.powi(2) - 1.) / k * hypot_dz;

    let glat = if deg { glat.to_degrees() } else { glat };
    let lon = if deg { lon.to_degrees() } else { lon };

    Some((glat, lon, alt))
}
