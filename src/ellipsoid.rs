#[derive(Clone, Debug, PartialEq)]
pub struct Ellipsoid {
    pub model: String,
    pub name: String,
    pub semimajor_axis: f64,
    pub semiminor_axis: f64,
    pub flattening: f64,
    pub thirdflattening: f64,
    pub eccentricity: f64,
}

impl Ellipsoid {
    pub fn new(semimajor_axis: f64, semiminor_axis: f64, name: &str, model: &str) -> Self {
        let flattening = (semimajor_axis - semiminor_axis) / semimajor_axis;
        assert!(flattening >= 0.0, "flattening must be >= 0");

        let thirdflattening = (semimajor_axis - semiminor_axis) / (semimajor_axis + semiminor_axis);
        let eccentricity = (2.0 * flattening - flattening.powi(2)).sqrt();

        Ellipsoid {
            model: model.to_string(),
            name: name.to_string(),
            semimajor_axis,
            semiminor_axis,
            flattening,
            thirdflattening,
            eccentricity,
        }
    }

    pub fn maupertuis() -> Ellipsoid {
        Ellipsoid::new(6397300.0, 6363806.283, "Maupertuis (1738)", "maupertuis")
    }
    pub fn plessis() -> Ellipsoid {
        Ellipsoid::new(6376523.0, 6355862.9333, "Plessis (1817)", "plessis")
    }
    pub fn everest1830() -> Ellipsoid {
        Ellipsoid::new(6377299.365, 6356098.359, "Everest (1830)", "everest1830")
    }
    pub fn everest1830m() -> Ellipsoid {
        Ellipsoid::new(
            6377304.063,
            6356103.039,
            "Everest 1830 Modified (1967)",
            "everest1830m",
        )
    }
    pub fn everest1967() -> Ellipsoid {
        Ellipsoid::new(
            6377298.556,
            6356097.55,
            "Everest 1830 (1967 Definition)",
            "everest1967",
        )
    }
    pub fn airy() -> Ellipsoid {
        Ellipsoid::new(6377563.396, 6356256.909, "Airy (1830)", "airy")
    }
    pub fn bessel() -> Ellipsoid {
        Ellipsoid::new(6377397.155, 6356078.963, "Bessel (1841)", "bessel")
    }
    pub fn clarke1866() -> Ellipsoid {
        Ellipsoid::new(6378206.4, 6356583.8, "Clarke (1866)", "clarke1866")
    }
    pub fn clarke1878() -> Ellipsoid {
        Ellipsoid::new(6378190.0, 6356456.0, "Clarke (1878)", "clarke1878")
    }
    pub fn clarke1860() -> Ellipsoid {
        Ellipsoid::new(6378249.145, 6356514.87, "Clarke (1880)", "clarke1860")
    }
    pub fn helmert() -> Ellipsoid {
        Ellipsoid::new(6378200.0, 6356818.17, "Helmert (1906)", "helmert")
    }
    pub fn hayford() -> Ellipsoid {
        Ellipsoid::new(6378388.0, 6356911.946, "Hayford (1910)", "hayford")
    }
    pub fn international1924() -> Ellipsoid {
        Ellipsoid::new(
            6378388.0,
            6356911.946,
            "International (1924)",
            "international1924",
        )
    }
    pub fn krassovsky1940() -> Ellipsoid {
        Ellipsoid::new(
            6378245.0,
            6356863.019,
            "Krassovsky (1940)",
            "krassovsky1940",
        )
    }
    pub fn wgs66() -> Ellipsoid {
        Ellipsoid::new(6378145.0, 6356759.769, "WGS66 (1966)", "wgs66")
    }
    pub fn australian() -> Ellipsoid {
        Ellipsoid::new(
            6378160.0,
            6356774.719,
            "Australian National (1966)",
            "australian",
        )
    }
    pub fn international1967() -> Ellipsoid {
        Ellipsoid::new(
            6378157.5,
            6356772.2,
            "New International (1967)",
            "international1967",
        )
    }
    pub fn grs67() -> Ellipsoid {
        Ellipsoid::new(6378160.0, 6356774.516, "GRS-67 (1967)", "grs67")
    }
    pub fn sa1969() -> Ellipsoid {
        Ellipsoid::new(6378160.0, 6356774.719, "South American (1969)", "sa1969")
    }
    pub fn wgs72() -> Ellipsoid {
        Ellipsoid::new(6378135.0, 6356750.52001609, "WGS-72 (1972)", "wgs72")
    }
    pub fn grs80() -> Ellipsoid {
        Ellipsoid::new(6378137.0, 6356752.31414036, "GRS-80 (1979)", "grs80")
    }
    pub fn wgs84() -> Ellipsoid {
        Ellipsoid::new(6378137.0, 6356752.31424518, "WGS-84 (1984)", "wgs84")
    }
    pub fn wgs84_mean() -> Ellipsoid {
        Ellipsoid::new(
            6371008.7714,
            6371008.7714,
            "WGS-84 (1984) Mean",
            "wgs84_mean",
        )
    }
    pub fn iers1989() -> Ellipsoid {
        Ellipsoid::new(6378136.0, 6356751.302, "IERS (1989)", "iers1989")
    }
    pub fn pz90_11() -> Ellipsoid {
        Ellipsoid::new(6378136.0, 6356751.3618, "ПЗ-90 (2011)", "pz90.11")
    }
    pub fn iers2003() -> Ellipsoid {
        Ellipsoid::new(6378136.6, 6356751.9, "IERS (2003)", "iers2003")
    }
    pub fn gsk2011() -> Ellipsoid {
        Ellipsoid::new(6378136.5, 6356751.758, "ГСК (2011)", "gsk2011")
    }
    pub fn mercury() -> Ellipsoid {
        Ellipsoid::new(2440500.0, 2438300.0, "Mercury", "mercury")
    }
    pub fn venus() -> Ellipsoid {
        Ellipsoid::new(6051800.0, 6051800.0, "Venus", "venus")
    }
    pub fn moon() -> Ellipsoid {
        Ellipsoid::new(1738100.0, 1736000.0, "Moon", "moon")
    }
    pub fn mars() -> Ellipsoid {
        Ellipsoid::new(3396900.0, 3376097.80585952, "Mars", "mars")
    }
    pub fn jupyter() -> Ellipsoid {
        Ellipsoid::new(71492000.0, 66770054.3475922, "Jupiter", "jupyter")
    }
    pub fn io() -> Ellipsoid {
        Ellipsoid::new(1829.7, 1815.8, "Io", "io")
    }
    pub fn saturn() -> Ellipsoid {
        Ellipsoid::new(60268000.0, 54364301.5271271, "Saturn", "saturn")
    }
    pub fn uranus() -> Ellipsoid {
        Ellipsoid::new(25559000.0, 24973000.0, "Uranus", "uranus")
    }
    pub fn neptune() -> Ellipsoid {
        Ellipsoid::new(24764000.0, 24341000.0, "Neptune", "neptune")
    }
    pub fn pluto() -> Ellipsoid {
        Ellipsoid::new(1188000.0, 1188000.0, "Pluto", "pluto")
    }

    pub fn get(name: &str) -> Option<Ellipsoid> {
        match name {
            "maupertuis" => Some(Ellipsoid::maupertuis()),
            "plessis" => Some(Ellipsoid::plessis()),
            "everest1830" => Some(Ellipsoid::everest1830()),
            "everest1830m" => Some(Ellipsoid::everest1830m()),
            "everest1967" => Some(Ellipsoid::everest1967()),
            "airy" => Some(Ellipsoid::airy()),
            "bessel" => Some(Ellipsoid::bessel()),
            "clarke1866" => Some(Ellipsoid::clarke1866()),
            "clarke1878" => Some(Ellipsoid::clarke1878()),
            "clarke1860" => Some(Ellipsoid::clarke1860()),
            "helmert" => Some(Ellipsoid::helmert()),
            "hayford" => Some(Ellipsoid::hayford()),
            "international1924" => Some(Ellipsoid::international1924()),
            "krassovsky1940" => Some(Ellipsoid::krassovsky1940()),
            "wgs66" => Some(Ellipsoid::wgs66()),
            "australian" => Some(Ellipsoid::australian()),
            "international1967" => Some(Ellipsoid::international1967()),
            "grs67" => Some(Ellipsoid::grs67()),
            "sa1969" => Some(Ellipsoid::sa1969()),
            "wgs72" => Some(Ellipsoid::wgs72()),
            "grs80" => Some(Ellipsoid::grs80()),
            "wgs84" => Some(Ellipsoid::wgs84()),
            "wgs84_mean" => Some(Ellipsoid::wgs84_mean()),
            "iers1989" => Some(Ellipsoid::iers1989()),
            "pz90.11" => Some(Ellipsoid::pz90_11()),
            "iers2003" => Some(Ellipsoid::iers2003()),
            "gsk2011" => Some(Ellipsoid::gsk2011()),
            "mercury" => Some(Ellipsoid::mercury()),
            "venus" => Some(Ellipsoid::venus()),
            "moon" => Some(Ellipsoid::moon()),
            "mars" => Some(Ellipsoid::mars()),
            "jupyter" => Some(Ellipsoid::jupyter()),
            "io" => Some(Ellipsoid::io()),
            "saturn" => Some(Ellipsoid::saturn()),
            "uranus" => Some(Ellipsoid::uranus()),
            "neptune" => Some(Ellipsoid::neptune()),
            "pluto" => Some(Ellipsoid::pluto()),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_ellipsoid() {
        use super::Ellipsoid;

        let ellipsoid = Ellipsoid::get("wgs84").unwrap();

        assert_eq!(ellipsoid.name, "WGS-84 (1984)");
        assert_eq!(ellipsoid.semimajor_axis, 6378137.0);
        assert_eq!(ellipsoid.semiminor_axis, 6356752.31424518);

        let ellipsoid = Ellipsoid::get("wgs84_mean").unwrap();

        assert_eq!(ellipsoid.name, "WGS-84 (1984) Mean");
        assert_eq!(ellipsoid.semimajor_axis, 6371008.7714);
        assert_eq!(ellipsoid.semiminor_axis, 6371008.7714);

        let ellipsoid1 = Ellipsoid::wgs72();
        let ellipsoid2 = Ellipsoid::get("wgs72").unwrap();

        assert_eq!(ellipsoid1, ellipsoid2);
    }
}
