use crate::{aer2enu, ecef2geodetic, enu2uvw, geodetic2ecef, Ellipsoid};
use std::f64::consts::FRAC_PI_2;

pub fn look_at_spheroid(
    lat0: f64,
    lon0: f64,
    h0: f64,
    az: f64,
    tilt: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let a = ell.semimajor_axis;
    let b = ell.semimajor_axis;
    let c = ell.semiminor_axis;

    let el = if deg { tilt - 90.0 } else { tilt - FRAC_PI_2 };

    let (e, n, u) = aer2enu(az, el, 1.0, deg);

    let (u, v, w) = enu2uvw(e, n, u, lat0, lon0, deg);
    let (x, y, z) = geodetic2ecef(lat0, lon0, h0, ell, deg)?;

    let value = -(a.powi(2)) * b.powi(2) * w * z
        - a.powi(2) * c.powi(2) * v * y
        - b.powi(2) * c.powi(2) * u * x;
    let radical = a.powi(2) * b.powi(2) * w.powi(2) + a.powi(2) * c.powi(2) * v.powi(2)
        - a.powi(2) * v.powi(2) * z.powi(2)
        + 2. * a.powi(2) * v * w * y * z
        - a.powi(2) * w.powi(2) * y.powi(2)
        + b.powi(2) * c.powi(2) * u.powi(2)
        - b.powi(2) * u.powi(2) * z.powi(2)
        + 2. * b.powi(2) * u * w * x * z
        - b.powi(2) * w.powi(2) * x.powi(2)
        - c.powi(2) * u.powi(2) * y.powi(2)
        + 2. * c.powi(2) * u * v * x * y
        - c.powi(2) * v.powi(2) * x.powi(2);

    let magnitude = a.powi(2) * b.powi(2) * w.powi(2)
        + a.powi(2) * c.powi(2) * v.powi(2)
        + b.powi(2) * c.powi(2) * u.powi(2);

    if radical < 0.0 {
        return None;
    }

    let d = (value - a * b * c * radical.sqrt()) / magnitude;

    if d < 0.0 {
        return None;
    }

    let (lat, lon, ..) = ecef2geodetic(x + d * u, y + d * v, z + d * w, ell, deg);

    Some((lat, lon, d))
}
