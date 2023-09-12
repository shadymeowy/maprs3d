use crate::{greenwichsrt, juliandate};

pub fn eci2ecef(x: f64, y: f64, z: f64, datetime: time::PrimitiveDateTime) -> (f64, f64, f64) {
    let gst = greenwichsrt(juliandate(datetime));

    rotate(x, y, z, gst)
}

pub fn ecef2eci(x: f64, y: f64, z: f64, datetime: time::PrimitiveDateTime) -> (f64, f64, f64) {
    let gst = greenwichsrt(juliandate(datetime));

    rotate(x, y, z, -gst)
}

fn rotate(x: f64, y: f64, z: f64, theta: f64) -> (f64, f64, f64) {
    let x_ = x * theta.cos() + y * theta.sin();
    let y_ = -x * theta.sin() + y * theta.cos();

    (x_, y_, z)
}
