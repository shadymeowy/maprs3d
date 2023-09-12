use crate::ellipsoid::Ellipsoid;
use crate::utils::is_close;
use crate::{ecef2eci, eci2ecef};
use std::f64::consts::FRAC_PI_2;

pub fn geodetic2ecef(
    lat: f64,
    lon: f64,
    alt: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (lat_, lon_) = if deg {
        (lat.to_radians(), lon.to_radians())
    } else {
        (lat, lon)
    };

    if lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let sin_lat = lat_.sin();
    let cos_lat = lat_.cos();
    let sin_lon = lon_.sin();
    let cos_lon = lon_.cos();

    let n = ell.semimajor_axis.powi(2)
        / (ell.semimajor_axis * cos_lat).hypot(ell.semiminor_axis * sin_lat);

    let x = (n + alt) * cos_lat * cos_lon;
    let y = (n + alt) * cos_lat * sin_lon;
    let z = (n * (ell.semiminor_axis / ell.semimajor_axis).powi(2) + alt) * sin_lat;

    Some((x, y, z))
}

pub fn ecef2geodetic(x: f64, y: f64, z: f64, ell: &Ellipsoid, deg: bool) -> (f64, f64, f64) {
    let r = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
    let e = (ell.semimajor_axis.powi(2) - ell.semiminor_axis.powi(2)).sqrt();

    let u = ((0.5 * (r.powi(2) - e.powi(2)))
        + (0.5 * ((r.powi(2) - e.powi(2)).hypot(2.0 * e * z))))
    .sqrt();
    let q = (x.powi(2) + y.powi(2)).sqrt();
    let hue = (u.powi(2) + e.powi(2)).sqrt();

    let beta = (hue / u * z / x.hypot(y)).atan();
    let beta = if beta.is_nan() {
        if is_close(z, 0., 1e-9, 0.) {
            0.
        } else if z > 0. {
            FRAC_PI_2
        } else {
            -FRAC_PI_2
        }
    } else {
        beta
    };

    let dbeta = ((ell.semiminor_axis * u - ell.semimajor_axis * hue + e.powi(2)) * beta.sin())
        / (ell.semimajor_axis * hue * (1.0 / beta.cos()) - e.powi(2) * beta.cos());

    let beta = beta + dbeta;

    let lat = (ell.semimajor_axis / ell.semiminor_axis * beta.tan()).atan();
    let lon = y.atan2(x);

    let alt =
        ((z - ell.semiminor_axis * beta.sin()).hypot(q - ell.semimajor_axis * beta.cos())).abs();

    let inside = (x.powi(2) / ell.semimajor_axis.powi(2)
        + y.powi(2) / ell.semimajor_axis.powi(2)
        + z.powi(2) / ell.semiminor_axis.powi(2))
        < 1.0;

    let alt = if inside { -alt } else { alt };

    if deg {
        let lat = lat.to_degrees();
        let lon = lon.to_degrees();
        (lat, lon, alt)
    } else {
        (lat, lon, alt)
    }
}

pub fn ecef2enuv(u: f64, v: f64, w: f64, lat0: f64, lon0: f64, deg: bool) -> (f64, f64, f64) {
    let (lat0_, lon0_) = if deg {
        (lat0.to_radians(), lon0.to_radians())
    } else {
        (lat0, lon0)
    };

    let sin_lat0 = lat0_.sin();
    let cos_lat0 = lat0_.cos();
    let sin_lon0 = lon0_.sin();
    let cos_lon0 = lon0_.cos();

    let t = cos_lon0 * u + sin_lon0 * v;
    let u = -sin_lon0 * u + cos_lon0 * v;
    let w = cos_lat0 * t + sin_lat0 * w;
    let v = -sin_lat0 * t + cos_lat0 * w;

    (u, v, w)
}

pub fn ecef2enu(
    x: f64,
    y: f64,
    z: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x0, y0, z0) = geodetic2ecef(lat0, lon0, h0, ell, deg)?;

    Some(uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg))
}

pub fn enu2uvw(east: f64, north: f64, up: f64, lat0: f64, lon0: f64, deg: bool) -> (f64, f64, f64) {
    let (lat0_, lon0_) = if deg {
        (lat0.to_radians(), lon0.to_radians())
    } else {
        (lat0, lon0)
    };

    let sin_lat0 = lat0_.sin();
    let cos_lat0 = lat0_.cos();
    let sin_lon0 = lon0_.sin();
    let cos_lon0 = lon0_.cos();

    let t = cos_lat0 * up - sin_lat0 * north;
    let w = sin_lat0 * up + cos_lat0 * north;
    let u = cos_lon0 * t - sin_lon0 * east;
    let v = sin_lon0 * t + cos_lon0 * east;

    (u, v, w)
}

pub fn uvw2enu(u: f64, v: f64, w: f64, lat0: f64, lon0: f64, deg: bool) -> (f64, f64, f64) {
    let (lat0_, lon0_) = if deg {
        (lat0.to_radians(), lon0.to_radians())
    } else {
        (lat0, lon0)
    };
    let sin_lat0 = lat0_.sin();
    let cos_lat0 = lat0_.cos();
    let sin_lon0 = lon0_.sin();
    let cos_lon0 = lon0_.cos();

    let t = cos_lon0 * u + sin_lon0 * v;
    let east = -sin_lon0 * u + cos_lon0 * v;
    let up = cos_lat0 * t + sin_lat0 * w;
    let north = -sin_lat0 * t + cos_lat0 * w;

    (east, north, up)
}

pub fn eci2geodetic(
    x: f64,
    y: f64,
    z: f64,
    t: time::PrimitiveDateTime,
    ell: &Ellipsoid,
    deg: bool,
) -> (f64, f64, f64) {
    let (xecef, yecef, zecef) = eci2ecef(x, y, z, t);

    ecef2geodetic(xecef, yecef, zecef, ell, deg)
}

pub fn geodetic2eci(
    lat: f64,
    lon: f64,
    alt: f64,
    t: time::PrimitiveDateTime,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x, y, z) = geodetic2ecef(lat, lon, alt, ell, deg)?;

    Some(ecef2eci(x, y, z, t))
}

pub fn enu2ecef(
    e1: f64,
    n1: f64,
    u1: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x0, y0, z0) = geodetic2ecef(lat0, lon0, h0, ell, deg)?;
    let (dx, dy, dz) = enu2uvw(e1, n1, u1, lat0, lon0, deg);
    Some((x0 + dx, y0 + dy, z0 + dz))
}
