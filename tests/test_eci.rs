use time_macros::datetime;

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
fn test_eci2ecef() {
    let eci = (-2981784., 5207055., 3161595.);
    let utc = datetime!(2019-1-4 12:00:00);
    let ecef = maprs3d::eci2ecef(eci.0, eci.1, eci.2, utc);

    let ecef2 = (-5.7627e6, -1.6827e6, 3.1560e6);

    assert!(is_close(ecef, ecef2, 0.025, 0.));
}

#[test]
fn test_ecef2eci() {
    let ecef = (-5762640., -1682738., 3156028.);
    let utc = datetime!(2019-1-4 12:00:00);
    let eci = maprs3d::ecef2eci(ecef.0, ecef.1, ecef.2, utc);

    let eci2 = (-2981810.6, 5207039.5, 3161595.1);
    assert!(is_close(eci, eci2, 0.01, 0.));
}

#[test]
fn test_eci2geodetic() {
    let ell = maprs3d::Ellipsoid::wgs84();
    let eci = (-2981784., 5207055., 3161595.);
    let utc = datetime!(2019-1-4 12:00:00);
    let geodetic = maprs3d::eci2geodetic(eci.0, eci.1, eci.2, utc, &ell, true);

    let geodetic2 = (27.880801, -163.722058, 408850.646);
    assert!(is_close(geodetic, geodetic2, 0.01, 0.));
}

#[test]
fn test_geodetic2eci() {
    let ell = maprs3d::Ellipsoid::wgs84();
    let geodetic = (27.880801, -163.722058, 408850.646);
    let utc = datetime!(2019-1-4 12:00:00);
    let eci = maprs3d::geodetic2eci(geodetic.0, geodetic.1, geodetic.2, utc, &ell, true).unwrap();

    let eci2 = (-2981784., 5207055., 3161595.);
    assert!(is_close(eci, eci2, 0.01, 0.));
}

#[test]
fn test_eci_aer() {
    let ell = maprs3d::Ellipsoid::wgs84();
    let eci = (4500000., -45000000., 3000000.);
    let t = datetime!(2022-1-2 03:04:05);
    let lla = (28., -80., 100.);
    let aer = maprs3d::eci2aer(eci.0, eci.1, eci.2, t, lla.0, lla.1, lla.2, &ell, true).unwrap();

    let aer2 = (314.9945, -53.0089, 5.026e7);
    assert!(is_close(aer, aer2, 0.01, 0.));

    let eci2 = maprs3d::aer2eci(aer.0, aer.1, aer.2, t, lla.0, lla.1, lla.2, &ell, true).unwrap();
    assert!(is_close(eci, eci2, 0.01, 0.));
}