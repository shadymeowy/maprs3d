use crate::{rcurve, vdist, Ellipsoid};

pub fn eqavol(ell: &Ellipsoid) -> f64 {
    let f = ell.flattening;

    ell.semimajor_axis * (1. - f / 3. - f.powi(2) / 9.)
}

pub fn authalic(ell: &Ellipsoid) -> f64 {
    let e = ell.eccentricity;

    if e > 0. {
        let f1 = ell.semimajor_axis.powi(2) / 2.;
        let f2 = (1. - e.powi(2)) / (2. * e);
        let f3 = 2. * e.atanh();

        (f1 * (1. + f2 * f3)).sqrt()
    } else {
        ell.semimajor_axis
    }
}

pub fn rectifying(ell: &Ellipsoid) -> f64 {
    ((ell.semimajor_axis.powf(3. / 2.) + ell.semiminor_axis.powf(3. / 2.)) / 2.).powf(2. / 3.)
}

pub fn euler(
    lat1: f64,
    lon1: f64,
    lat2: f64,
    lon2: f64,
    ell: &Ellipsoid,
    deg: bool,
) -> Option<f64> {
    let latmid = (lat1 + lat2) / 2.;
    let az = vdist(lat1, lon1, lat2, lon2, ell, deg)?.1;
    let az = if deg { az.to_radians() } else { az };

    let rho = rcurve::meridian(latmid, ell, deg)?;
    let nu = rcurve::transverse(latmid, ell, deg)?;

    let den = rho * az.sin().powi(2) + nu * az.cos().powi(2);

    Some((rho * nu) / den)
}

pub enum Method {
    Mean,
    Norm,
}

pub fn curve(lat: f64, ell: &Ellipsoid, deg: bool, method: Method) -> Option<f64> {
    let rho = rcurve::meridian(lat, ell, deg)?;
    let nu = rcurve::transverse(lat, ell, deg)?;

    match method {
        Method::Mean => Some((rho + nu) / 2.),
        Method::Norm => Some((rho * nu).sqrt()),
    }
}

pub fn triaxial(ell: &Ellipsoid, method: Method) -> f64 {
    match method {
        Method::Mean => (2. * ell.semimajor_axis + ell.semiminor_axis) / 3.,
        Method::Norm => (ell.semimajor_axis.powi(2) * ell.semiminor_axis).cbrt()
    }
}

pub fn biaxial(ell: &Ellipsoid, method: Method) -> f64 {
    match method {
        Method::Mean => (ell.semimajor_axis + ell.semiminor_axis) / 2.,
        Method::Norm => (ell.semimajor_axis * ell.semiminor_axis).sqrt()
    }
}