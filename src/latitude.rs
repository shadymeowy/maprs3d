use crate::{transverse, Ellipsoid};
use std::f64::{consts::FRAC_PI_2, INFINITY};

pub fn geoc2geod(
    geocentric_lat: f64,
    geocentric_distance: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<f64> {
    let geocentric_lat_ = if deg {
        geocentric_lat.to_radians()
    } else {
        geocentric_lat
    };

    if geocentric_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let r = geocentric_distance / ell.semimajor_axis;

    let geodetic_lat = geocentric_lat_
        + ((2. * geocentric_lat_).sin() / r) * ell.flattening
        + ((1. / r.powi(2) + 1. / (4. * r)) * (4. * geocentric_lat_).sin())
            * ell.flattening.powi(2);

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}

pub fn geodetic2geocentric(
    geodetic_lat: f64,
    alt_m: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<f64> {
    let geodetic_lat_ = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let r = transverse(geodetic_lat, ell, false)?;
    let geocentric_lat =
        ((1. - ell.eccentricity.powi(2) * (r / (r + alt_m))) * (geodetic_lat).tan()).atan();

    if deg {
        Some(geocentric_lat.to_degrees())
    } else {
        Some(geocentric_lat)
    }
}

pub use geodetic2geocentric as geod2geoc;

pub fn geocentric2geodetic(
    geocentric_lat: f64,
    alt_m: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<f64> {
    let geocentric_lat_ = if deg {
        geocentric_lat.to_radians()
    } else {
        geocentric_lat
    };

    if geocentric_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let r = transverse(geocentric_lat, ell, false)?;
    let geodetic_lat =
        (geocentric_lat.tan() / (1. - ell.eccentricity.powi(2) * (r / (r + alt_m)))).atan();

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}

pub fn geodetic2isometric(geodetic_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let geodetic_lat_ = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let e = ell.eccentricity;

    let mut isometric_lat = geodetic_lat.tan().asinh() - e * (e * geodetic_lat.sin()).atanh();
    let cos_lat = geodetic_lat.cos();

    if cos_lat.abs() <= 1e-9 {
        isometric_lat = INFINITY * geodetic_lat.signum();
    }

    if deg {
        Some(isometric_lat.to_degrees())
    } else {
        Some(isometric_lat)
    }
}

pub fn isometric2geodetic(isometric_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let isometric_lat_ = if deg {
        isometric_lat.to_radians()
    } else {
        isometric_lat
    };

    if isometric_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let conformal_lat = 2. * isometric_lat.exp().atan() - FRAC_PI_2;
    let geodetic_lat = conformal2geodetic(conformal_lat, ell, false)?;

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}

pub fn conformal2geodetic(conformal_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let conformal_lat_ = if deg {
        conformal_lat.to_radians()
    } else {
        conformal_lat
    };

    if conformal_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let e = ell.eccentricity;

    let f1 = e.powi(2) / 2. + 5. * e.powi(4) / 24. + e.powi(6) / 12. + 13. * e.powi(8) / 360.;
    let f2 = 7. * e.powi(4) / 48. + 2.9 * e.powi(6) / 240. + 811. * e.powi(8) / 11520.;
    let f3 = 7. * e.powi(6) / 120. + 81. * e.powi(8) / 1120.;
    let f4 = 4279. * e.powi(8) / 161280.;

    let geodetic_lat = conformal_lat_
        + f1 * (2. * conformal_lat_).sin()
        + f2 * (4. * conformal_lat_).sin()
        + f3 * (6. * conformal_lat_).sin()
        + f4 * (8. * conformal_lat_).sin();

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}

pub fn geodetic2conformal(geodetic_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let geodetic_lat_ = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let e = ell.eccentricity;

    let f1 = 1. - e * geodetic_lat.sin();
    let f2 = 1. + e * geodetic_lat.sin();
    let f3 = 1. - geodetic_lat.sin();
    let f4 = 1. + geodetic_lat.sin();

    let t = ((f4 / f3) * ((f1 / f2).powf(e))).sqrt();
    let conformal_lat = if t.is_nan() {
        FRAC_PI_2
    } else {
        2. * t.atan() - FRAC_PI_2
    };

    if deg {
        Some(conformal_lat.to_degrees())
    } else {
        Some(conformal_lat)
    }
}

pub fn geodetic2rectifying(geodetic_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let geodetic_lat_ = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let n = ell.thirdflattening;

    let f1 = 3. * n / 2. - 9. * n.powi(3) / 16.;
    let f2 = 15. * n.powi(2) / 16. - 15. * n.powi(4) / 32.;
    let f3 = 35. * n.powi(3) / 48.;
    let f4 = 315. * n.powi(4) / 512.;

    let rectifying_lat = geodetic_lat_ - f1 * (2. * geodetic_lat_).sin()
        + f2 * (4. * geodetic_lat_).sin()
        - f3 * (6. * geodetic_lat_).sin()
        + f4 * (8. * geodetic_lat_).sin();

    if deg {
        Some(rectifying_lat.to_degrees())
    } else {
        Some(rectifying_lat)
    }
}

pub fn rectifying2geodetic(rectifying_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let rectifying_lat_ = if deg {
        rectifying_lat.to_radians()
    } else {
        rectifying_lat
    };

    if rectifying_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let n = ell.thirdflattening;

    let f1 = 3. * n / 2. - 27. * n.powi(3) / 32.;
    let f2 = 21. * n.powi(2) / 16. - 55. * n.powi(4) / 32.;
    let f3 = 151. * n.powi(3) / 96.;
    let f4 = 1097. * n.powi(4) / 512.;

    let geodetic_lat = rectifying_lat_
        + f1 * (2. * rectifying_lat_).sin()
        + f2 * (4. * rectifying_lat_).sin()
        + f3 * (6. * rectifying_lat_).sin()
        + f4 * (8. * rectifying_lat_).sin();

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}

pub fn geodetic2authalic(geodetic_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let geodetic_lat_ = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let e = ell.eccentricity;

    let f1 = e.powi(2) / 3. + 31. * e.powi(4) / 180. + 59. * e.powi(6) / 560.;
    let f2 = 17. * e.powi(4) / 360. + 61. * e.powi(6) / 1260.;
    let f3 = 383. * e.powi(6) / 45360.;

    let authalic_lat = geodetic_lat - f1 * (2. * geodetic_lat).sin()
        + f2 * (4. * geodetic_lat).sin()
        - f3 * (6. * geodetic_lat).sin();

    if deg {
        Some(authalic_lat.to_degrees())
    } else {
        Some(authalic_lat)
    }
}

pub fn authalic2geodetic(
    authalic_lat: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<f64> {
    let authalic_lat_ = if deg {
        authalic_lat.to_radians()
    } else {
        authalic_lat
    };

    if authalic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let e = ell.eccentricity;

    let f1 = e.powi(2) / 3. + 31. * e.powi(4) / 180. + 517. * e.powi(6) / 5040.;
    let f2 = 23. * e.powi(4) / 360. + 251. * e.powi(6) / 3780.;
    let f3 = 761. * e.powi(6) / 45360.;

    let geodetic_lat = authalic_lat_
        + f1 * (2. * authalic_lat_).sin()
        + f2 * (4. * authalic_lat_).sin()
        + f3 * (6. * authalic_lat_).sin();

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}

pub fn geodetic2parametric(geodetic_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let geodetic_lat_ = if deg {
        geodetic_lat.to_radians()
    } else {
        geodetic_lat
    };

    if geodetic_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let parametric_lat = ((1. - ell.eccentricity.powi(2)).sqrt() * geodetic_lat.tan()).atan();

    if deg {
        Some(parametric_lat.to_degrees())
    } else {
        Some(parametric_lat)
    }
}

pub fn parametric2geodetic(parametric_lat: f64, ell: &Ellipsoid, deg: bool) -> Option<f64> {
    let parametric_lat_ = if deg {
        parametric_lat.to_radians()
    } else {
        parametric_lat
    };

    if parametric_lat_.abs() > FRAC_PI_2 {
        return None;
    }

    let geodetic_lat = (parametric_lat.tan() / (1. - ell.eccentricity.powi(2)).sqrt()).atan();

    if deg {
        Some(geodetic_lat.to_degrees())
    } else {
        Some(geodetic_lat)
    }
}
