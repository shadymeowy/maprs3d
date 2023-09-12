use time::PrimitiveDateTime;
use time_macros::datetime;

fn get_t0() -> PrimitiveDateTime {
    datetime!(2014-04-06 08:00:00)
}

#[test]
fn test_juliantime() {
    let jl = maprs3d::juliandate(get_t0());
    assert_eq!(jl, 2456753.8333333335);
}