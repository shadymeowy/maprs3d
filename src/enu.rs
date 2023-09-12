use std::f64::consts::TAU;

use crate::{ecef2geodetic, enu2ecef, geodetic2ecef, uvw2enu, Ellipsoid};

pub fn enu2aer(e: f64, n: f64, u: f64, deg: bool) -> (f64, f64, f64) {
    let e = if e.abs() < 1e-3 { 0. } else { e };
    let n = if n.abs() < 1e-3 { 0. } else { n };
    let u = if u.abs() < 1e-3 { 0. } else { u };

    let r = e.hypot(n);
    let slant_range = r.hypot(u);
    let elev = u.atan2(r);
    let mut az = n.atan2(e);

    if az < 0. {
        az = az.rem_euclid(TAU);
    }

    let az = if deg { az.to_degrees() } else { az };
    let elev = if deg { elev.to_degrees() } else { elev };
    
    (az, elev, slant_range)
}

pub fn aer2enu(az: f64, el: f64, srange: f64, deg: bool) -> (f64, f64, f64) {
    let az = if deg { az.to_radians() } else { az };
    let el = if deg { el.to_radians() } else { el };
    let srange = srange.abs();

    let r = srange * el.cos();
    let e = r * az.sin();
    let n = r * az.cos();
    let u = srange * el.sin();
    
    (e, n, u)
}

pub fn enu2geodetic(
    e: f64,
    n: f64,
    u: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x0, y0, z0) = enu2ecef(e, n, u, lat0, lon0, h0, ell, deg)?;

    Some(ecef2geodetic(x0, y0, z0, ell, deg))
}

pub fn geodetic2enu(
    lat: f64,
    lon: f64,
    h: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x1, y1, z1) = geodetic2ecef(lat, lon, h, ell, deg)?;
    let (x2, y2, z2) = geodetic2ecef(lat0, lon0, h0, ell, deg)?;

    Some(uvw2enu(x1 - x2, y1 - y2, z1 - z2, lat0, lon0, deg))
}
