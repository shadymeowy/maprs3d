use crate::{
    aer2enu, ecef2enu, ecef2geodetic, enu2aer, enu2uvw, geodetic2ecef, geodetic2enu, Ellipsoid, ecef2eci, eci2ecef,
};

pub fn ecef2aer(
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

    Some(enu2aer(e, n, u, deg))
}

pub fn geodetic2aer(
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

    Some(enu2aer(e, n, u, deg))
}

pub fn aer2geodetic(
    az: f64,
    el: f64,
    srange: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x, y, z) = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg)?;

    Some(ecef2geodetic(x, y, z, ell, deg))
}

pub fn eci2aer(
    x: f64,
    y: f64,
    z: f64,
    t: time::PrimitiveDateTime,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (xecef, yecef, zecef) = eci2ecef(x, y, z, t);

    ecef2aer(xecef, yecef, zecef, lat0, lon0, h0, ell, deg)
}

pub fn aer2eci(
    az: f64,
    el: f64,
    srange: f64,
    t: time::PrimitiveDateTime,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x, y, z) = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg)?;
    
    Some(ecef2eci(x, y, z, t))
}

pub fn aer2ecef(
    az: f64,
    el: f64,
    srange: f64,
    lat0: f64,
    lon0: f64,
    h0: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<(f64, f64, f64)> {
    let (x0, y0, z0) = geodetic2ecef(lat0, lon0, h0, ell, deg)?;
    let (e1, n1, u1) = aer2enu(az, el, srange, deg);
    let (dx, dy, dz) = enu2uvw(e1, n1, u1, lat0, lon0, deg);

    Some((x0 + dx, y0 + dy, z0 + dz))
}
