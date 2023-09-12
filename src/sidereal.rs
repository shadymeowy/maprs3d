use std::f64::consts::TAU;
use time;

pub fn datetime2sidereal(datetime: time::PrimitiveDateTime, lon_radians: f64) -> f64 {
    let jd = juliandate(datetime);
    let gst = greenwichsrt(jd);

    gst + lon_radians
}

pub fn juliandate(datetime: time::PrimitiveDateTime) -> f64 {
    let mut year = datetime.year() as f64;
    let mut month = datetime.month() as i32 as f64;
    let second = datetime.second() as f64;
    let minute = datetime.minute() as f64;
    let hour = datetime.hour() as f64;
    let day = datetime.day() as f64;

    if month < 3.0 {
        year -= 1.0;
        month += 12.0;
    }

    let a = (year / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    let c = ((second / 60.0 + minute) / 60.0 + hour) / 24.0;

    let jd = (365.25 * (year + 4716.0)).floor() + (30.6001 * (month + 1.0)).floor() + b + c
        - 1524.5
        + day;

    jd
}

pub fn greenwichsrt(jdate: f64) -> f64 {
    let tut1 = (jdate - 2451545.0) / 36525.0;

    let gmst_sec =
        67310.54841 + (876600. * 3600. + 8640184.812866) * tut1 + 0.093104 * tut1.powi(2)
            - 6.2e-6 * tut1.powi(3);

    (gmst_sec * TAU / 86400.0).rem_euclid(TAU)
}
