use crate::{
    aer2enu, ecef2enu, ecef2enuv, enu2aer, enu2ecef, enu2geodetic, geodetic2enu, Ellipsoid,
};

pub fn aer2ned(az: f64, el: f64, srange: f64, deg: bool) -> (f64, f64, f64) {
    let (e, n, u) = aer2enu(az, el, srange, deg);

    (n, e, -u)
}

pub fn ned2aer(n: f64, e: f64, d: f64, deg: bool) -> (f64, f64, f64) {
    enu2aer(e, n, -d, deg)
}

pub fn ned2geodetic(
    n: f64,
    e: f64,
    d: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    enu2geodetic(e, n, -d, lat0, lon0, h0, ell, deg)
}

pub fn ned2ecef(
    n: f64,
    e: f64,
    d: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (e, n, u) = enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg)?;

    Some((e, n, u))
}

pub fn ecef2ned(
    x: f64,
    y: f64,
    z: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (e, n, u) = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg)?;

    Some((n, e, -u))
}

pub fn geodetic2ned(
    lat: f64,
    lon: f64,
    h: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (e, n, u) = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg)?;

    Some((n, e, -u))
}

pub fn ecef2nedv(x: f64, y: f64, z: f64, lat0: f64, lon0: f64, deg: bool) -> (f64, f64, f64) {
    let (e, n, u) = ecef2enuv(x, y, z, lat0, lon0, deg);

    (n, e, -u)
}
