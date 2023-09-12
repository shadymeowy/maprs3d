pub fn anglesep(lon0: f64, lat0: f64, lon1: f64, lat1: f64, deg: bool) -> f64 {
    let lat0 = if deg { lat0.to_radians() } else { lat0 };
    let lon0 = if deg { lon0.to_radians() } else { lon0 };
    let lat1 = if deg { lat1.to_radians() } else { lat1 };

    let sep_rad = 2.
        * (haversine(lat0 - lat1) + lat0.cos() * lat1.cos() * haversine(lon0 - lon1))
            .sqrt()
            .asin();

    let sep_rad = if deg { sep_rad.to_degrees() } else { sep_rad };
    
    sep_rad
}

pub fn haversine(theta: f64) -> f64 {
    (1.0 - theta.cos()) / 2.0
}
