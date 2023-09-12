use time::{Date, Month, PrimitiveDateTime, Time};

fn get_t0() -> PrimitiveDateTime {
    PrimitiveDateTime::new(
        Date::from_calendar_date(2014, Month::April, 6).unwrap(),
        Time::from_hms(8, 0, 0).unwrap()
    )
}

#[test]
fn test_juliantime() {
    let jl = maprs3d::juliandate(get_t0());
    assert_eq!(jl, 2456753.8333333335);
}