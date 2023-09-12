#[test]
fn test_unit_vdist() {
    let params = vec![
        (0., 0., 0., 0., 0., 0.),
        (0., 0., 0., 90., 1.001875e7, 90.),
        (0., 0., 0., -90., 1.001875e7, 270.),
        (0., 0., 0., 180., 2.00375e7, 90.),
        (0., 0., 0., -180., 2.00375e7, 90.),
        (0., 0., 0., 4., 445277.96, 90.),
        (0., 0., 0., 5., 556597.45, 90.),
        (0., 0., 0., 6., 667916.94, 90.),
        (0., 0., 0., -6., 667916.94, 270.),
        (0., 0., 0., 7., 779236.44, 90.),
        (1e-16, 1e-16, 1e-16, 1., 111319.49, 90.),
        (90., 0., 0., 0., 1.00019657e7, 180.),
        (90., 0., -90., 0., 2.000393145e7, 180.),
    ];

    for (lat, lon, lat1, lon1, srange, az) in params {
        let (dist, az1) =
            maprs3d::vdist(lat, lon, lat1, lon1, &maprs3d::Ellipsoid::wgs84(), true).unwrap();
        dbg!(dist, srange);

        assert!(maprs3d::is_close(dist, srange, 0.005, 0.));
        assert!(maprs3d::is_close(az1, az, 1e-6, 1e-12));
    }
}

#[test]
fn test_identity() {
    let params = vec![(10., 20., 3e3, 38.), (0., 0., 0., 0.)];

    for (lat, lon, slantrange, az) in params {
        let (lat1, lon1) =
            maprs3d::vreckon(lat, lon, slantrange, az, &maprs3d::Ellipsoid::wgs84(), true).unwrap();

        let (dist, az1) =
            maprs3d::vdist(lat, lon, lat1, lon1, &maprs3d::Ellipsoid::wgs84(), true).unwrap();

        dbg!(dist, slantrange, az1, az);
        assert!(maprs3d::is_close(dist, slantrange, 1e-4, 1e-12)); // TODO
        assert!(maprs3d::is_close(az1, az, 1e-4, 1e-12)); // TODO
    }
}
