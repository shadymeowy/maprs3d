fn get_lla0() -> (f64, f64, f64) {
    (42., -82., 200.)
}

fn get_aer0() -> (f64, f64, f64) {
    (33., 70., 1000.)
}

fn is_close(abc: (f64, f64, f64), xyz: (f64, f64, f64), rel_tol: f64, abs_tol: f64) -> bool {
    // let rel_tol: f64 = 1e-9;
    // let abs_tol: f64 = 0.;

    let a = abc.0;
    let b = abc.1;
    let c = abc.2;
    let x = xyz.0;
    let y = xyz.1;
    let z = xyz.2;

    (a - x).abs() <= abs_tol.max(rel_tol * a.abs().max(x.abs()))
        && (b - y).abs() <= abs_tol.max(rel_tol * b.abs().max(y.abs()))
        && (c - z).abs() <= abs_tol.max(rel_tol * c.abs().max(z.abs()))
}

#[test]
fn test_ecef_ned() {
    let lla0 = get_lla0();
    let aer0 = get_aer0();
    let ell = maprs3d::Ellipsoid::wgs84();

    let enu = maprs3d::aer2enu(aer0.0, aer0.1, aer0.2, true);
    let ned = (enu.1, enu.0, -enu.2);
    let xyz = maprs3d::aer2ecef(aer0.0, aer0.1, aer0.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();

    let ned2 = maprs3d::ecef2ned(xyz.0, xyz.1, xyz.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    assert!(is_close(ned, ned2, 1e-9, 1e-12));
}

#[test]
fn test_ned_ecef() {
    let lla0 = get_lla0();
    let aer0 = get_aer0();
    let ell = maprs3d::Ellipsoid::wgs84();

    let enu = maprs3d::aer2enu(aer0.0, aer0.1, aer0.2, true);
    let ned = (enu.1, enu.0, -enu.2);
    let xyz = maprs3d::aer2ecef(aer0.0, aer0.1, aer0.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();

    let xyz2 = maprs3d::ned2ecef(ned.0, ned.1, ned.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    assert!(is_close(xyz, xyz2, 1e-9, 1e-12));
}


#[test]
fn test_enuv_nedv() {
    let lla0 = get_lla0();

    let (vx, vy, vz) = (5., 3., 2.);
    let (ve, vn, vu) = (5.368859646588048, 3.008520763668120, -0.352347711524077);
    
    let enuv = maprs3d::ecef2enuv(vx, vy, vz, lla0.0, lla0.1, true);
    assert!(is_close(enuv, (ve, vn, vu), 1e-9, 0.));

    let nedv = maprs3d::ecef2nedv(vx, vy, vz, lla0.0, lla0.1, true);
    assert!(is_close(nedv, (vn, ve, -vu), 1e-9, 0.));
}

#[test]
fn test_ned_geodetic() {
    let lla0 = get_lla0();
    let aer0 = get_aer0();
    let ell = maprs3d::Ellipsoid::wgs84();

    let lla1 = maprs3d::aer2geodetic(aer0.0, aer0.1, aer0.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    
    let enu3 = maprs3d::geodetic2enu(lla1.0, lla1.1, lla1.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    let ned3 = (enu3.1, enu3.0, -enu3.2);

    let ned2 = maprs3d::geodetic2ned(lla1.0, lla1.1, lla1.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    assert!(is_close(ned2, ned3, 1e-9, 1e-12));

    let lla2 = maprs3d::ned2geodetic(ned3.0, ned3.1, ned3.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    assert!(is_close(lla1, lla2, 1e-9, 1e-12));

    let lla2 = maprs3d::enu2geodetic(enu3.0, enu3.1, enu3.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    assert!(is_close(lla1, lla2, 1e-9, 1e-12));
}