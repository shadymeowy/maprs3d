use crate::Ellipsoid;
use std::f64::consts::FRAC_PI_2;

pub fn geocentric_radius(geodetic_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let geodetic_lat = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat.abs() > FRAC_PI_2 {
        return None;
    }

    let sin_lat = geodetic_lat.sin();
    let cos_lat = geodetic_lat.cos();

    let n = (ell.semimajor_axis.powi(2) * cos_lat).powi(2)
        + (ell.semiminor_axis.powi(2) * sin_lat).powi(2);
    let d = (ell.semimajor_axis * cos_lat).powi(2) + (ell.semiminor_axis * sin_lat).powi(2);

    Some((n / d).sqrt())
}

pub fn parallel(lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let lat = if deg { lat.to_radians() } else { lat };

    if lat.abs() > FRAC_PI_2 {
        return None;
    }

    Some(lat.cos() * transverse(lat, ell, false)?)
}

pub fn meridian(lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let lat = if deg { lat.to_radians() } else { lat };

    if lat.abs() > FRAC_PI_2 {
        return None;
    }

    let f1 = ell.semimajor_axis * (1. - ell.eccentricity.powi(2));
    let f2 = 1. - (ell.eccentricity * lat.sin()).powi(2);

    Some(f1 / (f2.powf(3. / 2.)))
}

pub fn transverse(lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let lat = if deg { lat.to_radians() } else { lat };

    if lat.abs() > FRAC_PI_2 {
        return None;
    }

    Some(ell.semimajor_axis / (1. - (ell.eccentricity * lat.sin()).powi(2)).sqrt())
}
