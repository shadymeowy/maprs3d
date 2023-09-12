#[test]
fn test_scalar_enu() {
    let ell = maprs3d::Ellipsoid::wgs84();
    let (lat0, lon0, h0) = (0., 90., -100.);
    let xyz = (0., ell.semimajor_axis, 50.);

    let enu = maprs3d::ecef2enu(xyz.0, xyz.1, xyz.2, lat0, lon0, h0, &ell, true).unwrap();
    let xyz2 = maprs3d::enu2ecef(enu.0, enu.1, enu.2, lat0, lon0, h0, &ell, true).unwrap();
    assert_eq!(xyz, xyz2);
}

#[test]
fn test_enu_ecef() {
    let ell = maprs3d::Ellipsoid::wgs84();

    let enu = (0., 0., 0.);
    let lla = (0., 0., 0.);
    let xyz = (ell.semimajor_axis, 0., 0.);

    let xyz2 = maprs3d::enu2ecef(enu.0, enu.1, enu.2, lla.0, lla.1, lla.2, &ell, true).unwrap();
    assert_eq!(xyz, xyz2);

    let rlla = (lla.0.to_radians(), lla.1.to_radians(), lla.2);
    let xyz3 = maprs3d::enu2ecef(enu.0, enu.1, enu.2, rlla.0, rlla.1, rlla.2, &ell, false).unwrap();
    assert_eq!(xyz, xyz3);

    let enu2 = maprs3d::ecef2enu(xyz.0, xyz.1, xyz.2, lla.0, lla.1, lla.2, &ell, true).unwrap();
    assert_eq!(enu, enu2);

    let enu2 = maprs3d::ecef2enu(xyz.0, xyz.1, xyz.2, rlla.0, rlla.1, rlla.2, &ell, false).unwrap();
    assert_eq!(enu, enu2);
}
