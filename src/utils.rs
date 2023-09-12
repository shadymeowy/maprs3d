pub fn cart2pol(x: f64, y: f64) -> (f64, f64) {
    let r = (x.powi(2) + y.powi(2)).sqrt();
    let theta = y.atan2(x);

    (r, theta)
}

pub fn pol2cart(theta: f64, rho: f64) -> (f64, f64) {
    let x = rho * theta.cos();
    let y = rho * theta.sin();
    
    (x, y)
}

pub fn cart2sph(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let hxy = x.hypot(y);
    let r = hxy.hypot(z);
    let el = z.atan2(hxy);
    let az = y.atan2(x);

    (az, el, r)
}

pub fn sph2cart(az: f64, el: f64, r: f64) -> (f64, f64, f64) {
    let rcos_theta = r * el.cos();
    let x = rcos_theta * az.cos();
    let y = rcos_theta * az.sin();
    let z = r * el.sin();

    (x, y, z)
}

pub fn sanitize(lat: f64, deg: bool) -> f64 {
    let lat = if deg { lat.to_radians() } else { lat };

    if lat.abs() > std::f64::consts::FRAC_PI_2 {
        panic!("-pi/2 <= latitude <= pi/2");
    }

    lat
}

pub fn is_close(a: f64, b: f64, rel_tol: f64, abs_tol: f64) -> bool {
    (a - b).abs() <= abs_tol.max(rel_tol * a.abs().max(b.abs()))
}
