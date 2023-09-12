use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq)]
pub struct Model {
    pub name: String,
    pub a: f64,
    pub b: f64,
}

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

    pub fn get_wgs84() -> Ellipsoid {
        Ellipsoid::new(6378137.0, 6356752.31424518, "WGS-84 (1984)", "wgs84")
    }
}

pub struct EllipsoidLibrary {
    models: HashMap<&'static str, Model>,
}

impl EllipsoidLibrary {
    pub fn new() -> Self {
        let models: HashMap<&'static str, Model> = [
            (
                "maupertuis",
                Model {
                    name: "Maupertuis (1738)".to_string(),
                    a: 6397300.0,
                    b: 6363806.283,
                },
            ),
            (
                "plessis",
                Model {
                    name: "Plessis (1817)".to_string(),
                    a: 6376523.0,
                    b: 6355862.9333,
                },
            ),
            (
                "everest1830",
                Model {
                    name: "Everest (1830)".to_string(),
                    a: 6377299.365,
                    b: 6356098.359,
                },
            ),
            (
                "everest1830m",
                Model {
                    name: "Everest 1830 Modified (1967)".to_string(),
                    a: 6377304.063,
                    b: 6356103.039,
                },
            ),
            (
                "everest1967",
                Model {
                    name: "Everest 1830 (1967 Definition)".to_string(),
                    a: 6377298.556,
                    b: 6356097.55,
                },
            ),
            (
                "airy",
                Model {
                    name: "Airy (1830)".to_string(),
                    a: 6377563.396,
                    b: 6356256.909,
                },
            ),
            (
                "bessel",
                Model {
                    name: "Bessel (1841)".to_string(),
                    a: 6377397.155,
                    b: 6356078.963,
                },
            ),
            (
                "clarke1866",
                Model {
                    name: "Clarke (1866)".to_string(),
                    a: 6378206.4,
                    b: 6356583.8,
                },
            ),
            (
                "clarke1878",
                Model {
                    name: "Clarke (1878)".to_string(),
                    a: 6378190.0,
                    b: 6356456.0,
                },
            ),
            (
                "clarke1860",
                Model {
                    name: "Clarke (1880)".to_string(),
                    a: 6378249.145,
                    b: 6356514.87,
                },
            ),
            (
                "helmert",
                Model {
                    name: "Helmert (1906)".to_string(),
                    a: 6378200.0,
                    b: 6356818.17,
                },
            ),
            (
                "hayford",
                Model {
                    name: "Hayford (1910)".to_string(),
                    a: 6378388.0,
                    b: 6356911.946,
                },
            ),
            (
                "international1924",
                Model {
                    name: "International (1924)".to_string(),
                    a: 6378388.0,
                    b: 6356911.946,
                },
            ),
            (
                "krassovsky1940",
                Model {
                    name: "Krassovsky (1940)".to_string(),
                    a: 6378245.0,
                    b: 6356863.019,
                },
            ),
            (
                "wgs66",
                Model {
                    name: "WGS66 (1966)".to_string(),
                    a: 6378145.0,
                    b: 6356759.769,
                },
            ),
            (
                "australian",
                Model {
                    name: "Australian National (1966)".to_string(),
                    a: 6378160.0,
                    b: 6356774.719,
                },
            ),
            (
                "international1967",
                Model {
                    name: "New International (1967)".to_string(),
                    a: 6378157.5,
                    b: 6356772.2,
                },
            ),
            (
                "grs67",
                Model {
                    name: "GRS-67 (1967)".to_string(),
                    a: 6378160.0,
                    b: 6356774.516,
                },
            ),
            (
                "sa1969",
                Model {
                    name: "South American (1969)".to_string(),
                    a: 6378160.0,
                    b: 6356774.719,
                },
            ),
            (
                "wgs72",
                Model {
                    name: "WGS-72 (1972)".to_string(),
                    a: 6378135.0,
                    b: 6356750.52001609,
                },
            ),
            (
                "grs80",
                Model {
                    name: "GRS-80 (1979)".to_string(),
                    a: 6378137.0,
                    b: 6356752.31414036,
                },
            ),
            (
                "wgs84",
                Model {
                    name: "WGS-84 (1984)".to_string(),
                    a: 6378137.0,
                    b: 6356752.31424518,
                },
            ),
            (
                "wgs84_mean",
                Model {
                    name: "WGS-84 (1984) Mean".to_string(),
                    a: 6371008.7714,
                    b: 6371008.7714,
                },
            ),
            (
                "iers1989",
                Model {
                    name: "IERS (1989)".to_string(),
                    a: 6378136.0,
                    b: 6356751.302,
                },
            ),
            (
                "pz90.11",
                Model {
                    name: "ПЗ-90 (2011)".to_string(),
                    a: 6378136.0,
                    b: 6356751.3618,
                },
            ),
            (
                "iers2003",
                Model {
                    name: "IERS (2003)".to_string(),
                    a: 6378136.6,
                    b: 6356751.9,
                },
            ),
            (
                "gsk2011",
                Model {
                    name: "ГСК (2011)".to_string(),
                    a: 6378136.5,
                    b: 6356751.758,
                },
            ),
            (
                "mercury",
                Model {
                    name: "Mercury".to_string(),
                    a: 2440500.0,
                    b: 2438300.0,
                },
            ),
            (
                "venus",
                Model {
                    name: "Venus".to_string(),
                    a: 6051800.0,
                    b: 6051800.0,
                },
            ),
            (
                "moon",
                Model {
                    name: "Moon".to_string(),
                    a: 1738100.0,
                    b: 1736000.0,
                },
            ),
            (
                "mars",
                Model {
                    name: "Mars".to_string(),
                    a: 3396900.0,
                    b: 3376097.80585952,
                },
            ),
            (
                "jupyter",
                Model {
                    name: "Jupiter".to_string(),
                    a: 71492000.0,
                    b: 66770054.3475922,
                },
            ),
            (
                "io",
                Model {
                    name: "Io".to_string(),
                    a: 1829.7,
                    b: 1815.8,
                },
            ),
            (
                "saturn",
                Model {
                    name: "Saturn".to_string(),
                    a: 60268000.0,
                    b: 54364301.5271271,
                },
            ),
            (
                "uranus",
                Model {
                    name: "Uranus".to_string(),
                    a: 25559000.0,
                    b: 24973000.0,
                },
            ),
            (
                "neptune",
                Model {
                    name: "Neptune".to_string(),
                    a: 24764000.0,
                    b: 24341000.0,
                },
            ),
            (
                "pluto",
                Model {
                    name: "Pluto".to_string(),
                    a: 1188000.0,
                    b: 1188000.0,
                },
            ),
        ]
        .iter()
        .cloned()
        .collect();

        EllipsoidLibrary { models }
    }

    pub fn get_ellipsoid(&self, name: &str) -> Option<Ellipsoid> {
        if let Some(model) = self.models.get(name) {
            Some(Ellipsoid::new(model.a, model.b, &model.name, name))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_ellipsoid() {
        use super::EllipsoidLibrary;

        let library = EllipsoidLibrary::new();

        let ellipsoid = library.get_ellipsoid("wgs84").unwrap();

        assert_eq!(ellipsoid.name, "WGS-84 (1984)");
        assert_eq!(ellipsoid.semimajor_axis, 6378137.0);
        assert_eq!(ellipsoid.semiminor_axis, 6356752.31424518);

        let ellipsoid = library.get_ellipsoid("wgs84_mean").unwrap();

        assert_eq!(ellipsoid.name, "WGS-84 (1984) Mean");
        assert_eq!(ellipsoid.semimajor_axis, 6371008.7714);
        assert_eq!(ellipsoid.semiminor_axis, 6371008.7714);
    }
}
