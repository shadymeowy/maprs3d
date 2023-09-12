use std::f64::consts::FRAC_2_PI;

use crate::datetime2sidereal;

pub fn azel2radec(
    az: f64,
    el: f64,
    lat: f64,
    lon: f64,
    deg: bool,
    datetime: time::PrimitiveDateTime,
) -> Option<(f64, f64)> {
    let (az_, el_, lat_, lon_) = if deg {
        (
            az.to_radians(),
            el.to_radians(),
            lat.to_radians(),
            lon.to_radians(),
        )
    } else {
        (az, el, lat, lon)
    };

    if lat_.abs() > FRAC_2_PI {
        return None;
    }

    let dec = el_.sin() * lat_.sin() + el_.cos() * lat_.cos() * az_.sin();

    let lha = (-(az_.sin() * el_.cos()) / dec.cos())
        .atan2((el_.sin() - lat_.sin() * dec.sin()) / (dec.cos() * lat_.cos()));

    let lst = datetime2sidereal(datetime, lon_);
    
    Some(((lst-lha).to_degrees().rem_euclid(360.0), dec.to_degrees()))
}

pub fn radec2azel(
    ra: f64,
    dec: f64,
    lat: f64,
    lon: f64,
    deg: bool,
    datetime: time::PrimitiveDateTime,
) -> Option<(f64, f64)> {
    let (ra_, dec_, lat_, lon_) = if deg {
        (
            ra.to_radians(),
            dec.to_radians(),
            lat.to_radians(),
            lon.to_radians(),
        )
    } else {
        (ra, dec, lat, lon)
    };

    if lat_.abs() > FRAC_2_PI {
        return None;
    }

    let ha = (datetime2sidereal(datetime, lon_) - ra_).to_radians();

    let el = dec_.sin() * lat_.sin() + dec_.cos() * lat_.cos() * ha.sin();

    let az = (-(ha.sin() * dec_.cos()) / el.cos())
        .atan2((dec_.sin() - lat_.sin() * el.sin()) / (el.cos() * lat.cos()));

    Some((az.to_degrees().rem_euclid(360.0), el.to_degrees()))
}