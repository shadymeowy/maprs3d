use std::f64::consts::{FRAC_2_PI, TAU};

use crate::datetime2sidereal;

pub fn azel2radec(
    az: f64,
    el: f64,
    lat: f64,
    lon: f64,
    deg: bool,
    datetime: time::PrimitiveDateTime,
) -> Option<(f64, f64)> {
    let az = if deg { az.to_radians() } else { az };
    let el = if deg { el.to_radians() } else { el };
    let lat = if deg { lat.to_radians() } else { lat };
    let lon = if deg { lon.to_radians() } else { lon };

    if lat.abs() > FRAC_2_PI {
        return None;
    }

    let dec = el.sin() * lat.sin() + el.cos() * lat.cos() * az.sin();

    let lha = (-(az.sin() * el.cos()) / dec.cos())
        .atan2((el.sin() - lat.sin() * dec.sin()) / (dec.cos() * lat.cos()));

    let lst = datetime2sidereal(datetime, lon);
    let ra = (lst - lha).rem_euclid(TAU);

    let ra = if deg { ra.to_degrees() } else { ra };
    let dec = if deg { dec.to_degrees() } else { dec };

    Some((ra, dec))
}

pub fn radec2azel(
    ra: f64,
    dec: f64,
    lat: f64,
    lon: f64,
    deg: bool,
    datetime: time::PrimitiveDateTime,
) -> Option<(f64, f64)> {
    let ra = if deg { ra.to_radians() } else { ra };
    let dec = if deg { dec.to_radians() } else { dec };
    let lat = if deg { lat.to_radians() } else { lat };
    let lon = if deg { lon.to_radians() } else { lon };

    if lat.abs() > FRAC_2_PI {
        return None;
    }

    let ha = (datetime2sidereal(datetime, lon) - ra).to_radians();

    let el = dec.sin() * lat.sin() + dec.cos() * lat.cos() * ha.sin();

    let az = (-(ha.sin() * dec.cos()) / el.cos())
        .atan2((dec.sin() - lat.sin() * el.sin()) / (el.cos() * lat.cos()))
        .rem_euclid(TAU);

    let az = if deg { az.to_degrees() } else { az };
    let el = if deg { el.to_degrees() } else { el };

    Some((az, el))
}
