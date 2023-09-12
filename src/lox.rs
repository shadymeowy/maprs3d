use std::f64::consts::PI;

use crate::{
    authalic2geodetic, geodetic2authalic, geodetic2isometric, geodetic2rectifying, rcurve,
    rectifying2geodetic, rsphere, sph2cart, Ellipsoid, cart2sph,
};

pub fn meridian_dist(lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    meridian_arc(0., lat, ell, deg)
}

pub fn meridian_arc(lat1: f64, lat2: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let lat1 = if deg { lat1.to_radians() } else { lat1 };
    let lat2 = if deg { lat2.to_radians() } else { lat2 };

    let rlat1 = geodetic2rectifying(lat1, ell, false)?;
    let rlat2 = geodetic2rectifying(lat2, ell, false)?;

    Some(rsphere::rectifying(ell) * (rlat2 - rlat1))
}

pub fn loxodrome_inverse(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64)> {
    let lat1 = if deg { lat1.to_radians() } else { lat1 };
    let lat2 = if deg { lat2.to_radians() } else { lat2 };
    let lon1 = if deg { lon1.to_radians() } else { lon1 };
    let lon2 = if deg { lon2.to_radians() } else { lon2 };

    let lati1 = geodetic2isometric(lat1, ell, false)?;
    let lati2 = geodetic2isometric(lat2, ell, false)?;

    let disolat = lati2 - lati1;
    let dlon = lon2 - lon1;

    let mut az12 = dlon.atan2(disolat);
    let aux = az12.cos().abs();

    let dist = if aux < 1e-9 {
        departure(lon2, lon1, lat1, ell, false)?
    } else {
        meridian_arc(lat1, lat2, ell, false)? / aux
    };

    if deg {
        az12 = az12.to_degrees().rem_euclid(360.);
    }

    Some((dist, az12))
}

pub fn loxodrome_direct(
    lat1: f64,
    lon1: f64,
    rng: f64,
    az12: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64)> {
    let lat1 = if deg { lat1.to_radians() } else { lat1 };
    let lon1 = if deg { lon1.to_radians() } else { lon1 };
    let az12 = if deg { az12.to_radians() } else { az12 };

    if rng < 0. || lat1.abs() > PI / 2. {
        return None;
    }

    let reclat = geodetic2rectifying(lat1, ell, false)?;

    let coaz = az12.cos();
    let rlat2 = reclat + (rng / rsphere::rectifying(ell)) * coaz;
    let lat2 = rectifying2geodetic(rlat2, ell, false)?;

    let newiso = geodetic2isometric(lat2, ell, false)?;
    let iso = geodetic2isometric(lat1, ell, false)?;

    let dlon = if coaz.abs() < 1e-9 {
        az12.tan() * (newiso - iso)
    } else {
        (PI - az12).signum() * rng / rcurve::parallel(lat1, ell, false)?
    };

    let lon2 = lon1 + dlon;

    if deg {
        Some((lat2.to_degrees(), lon2.to_degrees()))
    } else {
        Some((lat2, lon2))
    }
}

pub fn departure(lon1: f64, lon2: f64, lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let lon1 = if deg { lon1.to_radians() } else { lon1 };
    let lon2 = if deg { lon2.to_radians() } else { lon2 };
    let lat = if deg { lat.to_radians() } else { lat };

    Some(rcurve::parallel(lat, ell, false)? * (lon2 - lon1).abs().rem_euclid(PI))
}

pub fn meanm(lats: &[f64], lons: &[f64], ell: &Ellipsoid, deg: bool) -> Option<(f64, f64)> {
    let lats: Vec<f64> = if deg {
        lats.iter().map(|&x| x.to_radians()).collect()
    } else {
        lats.to_vec()
    };
    let lons = if deg {
        lons.iter().map(|&x| x.to_radians()).collect()
    } else {
        lons.to_vec()
    };

    let latas: Option<Vec<_>> = lats
        .iter()
        .map(|&x| geodetic2authalic(x, ell, false))
        .collect();
    let latas = latas?;

    let xyz: Vec<_> = latas
        .iter()
        .zip(lons.iter())
        .map(|(&lat, &lon)| sph2cart(lon, lat, 1.))
        .collect();
    
    let x = xyz.iter().map(|&x| x.0).sum();
    let y = xyz.iter().map(|&x| x.1).sum();
    let z = xyz.iter().map(|&x| x.2).sum();

    let (lonbar, latbar, ..) = cart2sph(x, y, z);
    let latbar = authalic2geodetic(latbar, ell, false)?;

    if deg {
        Some((latbar.to_degrees(), lonbar.to_degrees()))
    } else {
        Some((latbar, lonbar))
    }
}
