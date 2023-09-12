#[test]
fn test_track2() {
    let ell = maprs3d::Ellipsoid::wgs84();
    let llas = maprs3d::track2(40., 80., 65., -148., &ell, true, 3).unwrap();

    let llas2 = vec![
        (40., 80.),
        (69.633139886, 113.06849104),
        (65., -148.)
    ];
    
    for (lla, lla2) in llas.iter().zip(llas2.iter()) {
        assert!(maprs3d::is_close(lla.0, lla2.0, 1e-6, 1e-12));
        assert!(maprs3d::is_close(lla.1, lla2.1, 1e-6, 1e-12));
    }
}