pub fn anglesep(lon0: f64, lat0: f64, lon1: f64, lat1: f64, deg: bool) -> f64 {
    let (lon0_, lat0_, lon1_, lat1_) = {
        if deg {
            (
                lon0.to_radians(),
                lat0.to_radians(),
                lon1.to_radians(),
                lat1.to_radians(),
            )
        } else {
            (lon0, lat0, lon1, lat1)
        }
    };

    let sep_rad = 2.
        * (haversine(lat0_ - lat1_) + lat0.cos() * lat1.cos() * haversine(lon0_ - lon1_))
            .sqrt()
            .asin();

    if deg {
        sep_rad.to_degrees()
    } else {
        sep_rad
    }
}

pub fn haversine(theta: f64) -> f64 {
    (1.0 - theta.cos()) / 2.0
}
