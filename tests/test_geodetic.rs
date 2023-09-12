use std::{vec, f64::NAN};

fn get_lla0() -> (f64, f64, f64) {
    (42., -82., 200.)
}

fn get_rlla0() -> (f64, f64, f64) {
    let lla0 = get_lla0();
    (lla0.0.to_radians(), lla0.1.to_radians(), lla0.2)
}

fn get_xyz0() -> (f64, f64, f64) {
    (660675.2518247, -4700948.68316, 4245737.66222)
}

fn get_xyzlla() -> Vec<((f64, f64, f64), (f64, f64, f64))> {
    let ell = maprs3d::Ellipsoid::wgs84();
    let a = ell.semimajor_axis;
    let b = ell.semiminor_axis;
    vec![
        ((a, 0., 0.), (0., 0., 0.)),
        ((a - 1., 0., 0.), (0., 0., -1.)),
        ((a + 1., 0., 0.), (0., 0., 1.)),
        ((0.1 * a, 0., 0.), (0., 0., -0.9 * a)),
        ((0.001 * a, 0., 0.), (0., 0., -0.999 * a)),
        ((0., a, 0.), (0., 90., 0.)),
        ((0., a - 1., 0.), (0., 90., -1.)),
        ((0., a + 1., 0.), (0., 90., 1.)),
        ((0., 0.1 * a, 0.), (0., 90., -0.9 * a)),
        ((0., 0.001 * a, 0.), (0., 90., -0.999 * a)),
        ((0., 0., b), (90., 0., 0.)),
        ((0., 0., b + 1.), (90., 0., 1.)),
        ((0., 0., b - 1.), (90., 0., -1.)),
        ((0., 0., 0.1 * b), (90., 0., -0.9 * b)),
        ((0., 0., 0.001 * b), (90., 0., -0.999 * b)),
        ((0., 0., b - 1.), (89.999999, 0., -1.)),
        ((0., 0., b - 1.), (89.99999, 0., -1.)),
        ((0., 0., -b + 1.), (-90., 0., -1.)),
        ((0., 0., -b + 1.), (-89.999999, 0., -1.)),
        ((0., 0., -b + 1.), (-89.99999, 0., -1.)),
        ((-a + 1., 0., 0.), (0., 180., -1.)),
    ]
}

fn get_llaxyz() -> Vec<((f64, f64, f64), (f64, f64, f64))> {
    let ell = maprs3d::Ellipsoid::wgs84();
    let a = ell.semimajor_axis;
    let b = ell.semiminor_axis;

    vec![
        ((0., 0., -1.), (a - 1., 0., 0.)),
        ((0., 90., -1.), (0., a - 1., 0.)),
        ((0., -90., -1.), (0., -a + 1., 0.)),
        ((90., 0., -1.), (0., 0., b - 1.)),
        ((90., 15., -1.), (0., 0., b - 1.)),
        ((-90., 0., -1.), (0., 0., -b + 1.)),
    ]
}

fn get_aerllalla0() -> Vec<((f64, f64, f64), (f64, f64, f64), (f64, f64, f64))> {
    vec![
        ((33., 77., 1000.), (42.0016981935, -81.99852, 1174.374035), (42., -82., 200.)),
        ((0., 90., 10000.), (0., 0., 10000.), (0., 0., 0.))
    ]
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
fn test_scalar_geodetic2ecef() {
    let ell = maprs3d::Ellipsoid::wgs84();

    let lla = get_lla0();
    let xyz = maprs3d::geodetic2ecef(lla.0, lla.1, lla.2, &ell, true).unwrap();
    let lla1 = maprs3d::ecef2geodetic(xyz.0, xyz.1, xyz.2, &ell, true);
    assert!(is_close(lla, lla1, 1e-4, 0.));
}

#[test]
fn test_scalar_ecef2geodetic() {
    let ell = maprs3d::Ellipsoid::wgs84();

    let xyz = get_xyz0();
    let lla = maprs3d::ecef2geodetic(xyz.0, xyz.1, xyz.2, &ell, true);
    let xyz1 = maprs3d::geodetic2ecef(lla.0, lla.1, lla.2, &ell, true).unwrap();
    assert!(is_close(xyz, xyz1, 1e-4, 0.));
}

#[test]
fn test_ecef() {
    let ell = maprs3d::Ellipsoid::wgs84();
    let a = ell.semimajor_axis;

    let lla = get_lla0();
    let rlla = get_rlla0();
    let xyz1 = maprs3d::geodetic2ecef(lla.0, lla.1, lla.2, &ell, true).unwrap();
    let xyz2 = maprs3d::geodetic2ecef(rlla.0, rlla.1, rlla.2, &ell, false).unwrap();
    assert_eq!(xyz1, xyz2);

    assert!(maprs3d::geodetic2ecef(91., 0., 0., &ell, true).is_none());

    let lla2 = maprs3d::ecef2geodetic(xyz1.0, xyz1.1, xyz1.2, &ell, true);
    let lla3 = maprs3d::ecef2geodetic(xyz2.0, xyz2.1, xyz2.2, &ell, true);
    assert!(is_close(lla, lla2, 1e-9, 0.));
    assert!(is_close(lla, lla3, 1e-9, 0.));

    let lla = (0., 45., -1.);
    let s2 = 2_f64.sqrt();
    let xyz = maprs3d::ecef2geodetic((a - 1.) / s2, (a - 1.) / s2, 0., &ell, true);
    assert!(is_close(lla, xyz, 1e-9, 0.));
}

#[test]
fn test_geodetic2ecef() {
    let ell = maprs3d::Ellipsoid::wgs84();

    let llaxyz = get_llaxyz();
    for (lla, xyz) in llaxyz {
        let xyz1 = maprs3d::geodetic2ecef(lla.0, lla.1, lla.2, &ell, true).unwrap();
        assert!(is_close(xyz, xyz1, 1e-9, 1e-6));
    }
}

#[test]
fn test_ecef2geodetic() {
    let ell = maprs3d::Ellipsoid::wgs84();

    let xyzlla = get_xyzlla();
    for (xyz, lla) in xyzlla {
        let lla1 = maprs3d::ecef2geodetic(xyz.0, xyz.1, xyz.2, &ell, true);
        assert!(is_close(lla, lla1, 1e-6, 1e-9));
    }
}

#[test]
fn test_aer_geodetic() {
    let ell = maprs3d::Ellipsoid::wgs84();

    for (aer, lla, lla0) in get_aerllalla0() {
        let lla2 =
            maprs3d::aer2geodetic(aer.0, aer.1, aer.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();

        assert!(is_close(lla, lla2, 1e-6, 1e-12));

        let raer = (aer.0.to_radians(), aer.1.to_radians(), aer.2);
        let rlla = (lla.0.to_radians(), lla.1.to_radians(), lla.2);
        let rlla0 = (lla0.0.to_radians(), lla0.1.to_radians(), lla0.2);

        let rlla2 = maprs3d::aer2geodetic(raer.0, raer.1, raer.2, rlla0.0, rlla0.1, rlla0.2, &ell, false)
            .unwrap();

        assert!(is_close(rlla, rlla2, 1e-6, 1e-12));

        let aer2 = maprs3d::geodetic2aer(lla.0, lla.1, lla.2, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
        assert!(is_close(aer, aer2, 1e-3, 1e-12));

        let raer2 = maprs3d::geodetic2aer(rlla.0, rlla.1, rlla.2, rlla0.0, rlla0.1, rlla0.2, &ell, false)
            .unwrap();
        assert!(is_close(raer, raer2, 1e-3, 1e-12));
    }
}

#[test]
fn test_scalar_nan() {
    let lla0 = get_lla0();
    let ell = maprs3d::Ellipsoid::wgs84();

    let aer = maprs3d::geodetic2aer(NAN, NAN, NAN, lla0.0, lla0.1, lla0.2, &ell, true).unwrap();
    assert!(aer.0.is_nan());
    assert!(aer.1.is_nan());
    assert!(aer.2.is_nan());
}
