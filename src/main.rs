extern crate simd;

use simd::x86::sse2::f64x2;
use std::env;
use std::f64::consts::PI;

const BODIES: usize = 5;
const SOLAR_MASS: f64 = 4. * PI * PI;
const DAYS_PER_YEAR: f64 = 365.24;
const TIME_STEP: f64 = 0.01;

fn main() {
    let steps: u32 = env::args()
        .nth(1)
        .and_then(|a| a.parse().ok())
        .unwrap_or(1000);
    let mut bodies = Bodies::new();
    println!("{:.9}", bodies.get_energy());
    bodies.step_time(steps);
    println!("{:.9}", bodies.get_energy());
}

struct Bodies {
    bodies: [Body; BODIES]
}

impl Bodies {
    fn new() -> Self {
        let mut bodies = [
            Body::new_sun_no_velocity(),
            Body::new_jupiter(),
            Body::new_saturn(),
            Body::new_uranus(),
            Body::new_neptune()
        ];
        for i in 1..BODIES {
            bodies[0].vxy = bodies[0].vxy - bodies[i].vxy * f64x2::splat(bodies[i].mass / SOLAR_MASS);
            bodies[0].vz -= bodies[i].vz * bodies[i].mass / SOLAR_MASS;
        }
        Bodies { bodies }
    }

    fn get_energy(&self) -> f64 {
        let mut e = 0.;
        for i in 0..BODIES {
            let body_i = &self.bodies[i];
            let vxy_pow = body_i.vxy * body_i.vxy;
            e += 0.5 * body_i.mass * (
                vxy_pow.extract(0) +
                vxy_pow.extract(1) +
                body_i.vz.powi(2));
            for j in i + 1..BODIES {
                let body_j = &self.bodies[j];
                let mut dxy = body_i.xy - body_j.xy;
                dxy = dxy * dxy;
                let d =
                    dxy.extract(0) +
                    dxy.extract(1) +
                    (body_i.z - body_j.z).powi(2);
                e -= body_i.mass * body_j.mass / d.sqrt();
            }
        }
        e
    }

    fn step_time(&mut self, steps: u32) {
        for _ in 0..steps {
            for i in 0..BODIES {
                for j in i + 1..BODIES {
                    let dxy = self.bodies[i].xy - self.bodies[j].xy;
                    let dz = self.bodies[i].z - self.bodies[j].z;

                    let dxy_2 = dxy * dxy;
                    let d = dxy_2.extract(0) + dxy_2.extract(1) + dz * dz;
                    let mag = TIME_STEP / (d.sqrt() * d);

                    let m_j_mag = self.bodies[j].mass * mag;
                    self.bodies[i].vxy = self.bodies[i].vxy - dxy * f64x2::splat(m_j_mag);
                    self.bodies[i].vz -= dz * m_j_mag;

                    let m_i_mag = self.bodies[i].mass * mag;
                    self.bodies[j].vxy = self.bodies[j].vxy + dxy * f64x2::splat(m_i_mag);
                    self.bodies[j].vz += dz * m_i_mag;
                }
                self.bodies[i].xy = self.bodies[i].xy + self.bodies[i].vxy * f64x2::splat(TIME_STEP);
                self.bodies[i].z += self.bodies[i].vz * TIME_STEP;
            }
        }
    }
}

struct Body {
    xy: f64x2,
    z: f64,
    vxy: f64x2,
    vz: f64,
    mass: f64,
}

impl Body {
    fn new_jupiter() -> Self {
        Body {
            xy: f64x2::new(4.84143144246472090e+00, -1.16032004402742839e+00),
            z: -1.03622044471123109e-01,
            vxy: f64x2::new(1.66007664274403694e-03 * DAYS_PER_YEAR, 7.69901118419740425e-03 * DAYS_PER_YEAR),
            vz: -6.90460016972063023e-05 * DAYS_PER_YEAR,
            mass: 9.54791938424326609e-04 * SOLAR_MASS,
        }
    }

    fn new_saturn() -> Self {
        Body {
            xy: f64x2::new(8.34336671824457987e+00, 4.12479856412430479e+00),
            z: -4.03523417114321381e-01,
            vxy: f64x2::new(-2.76742510726862411e-03 * DAYS_PER_YEAR, 4.99852801234917238e-03 * DAYS_PER_YEAR),
            vz: 2.30417297573763929e-05 * DAYS_PER_YEAR,
            mass: 2.85885980666130812e-04 * SOLAR_MASS,
        }
    }

    fn new_uranus() -> Self {
        Body {
            xy: f64x2::new(1.28943695621391310e+01, -1.51111514016986312e+01),
            z: -2.23307578892655734e-01,
            vxy: f64x2::new(2.96460137564761618e-03 * DAYS_PER_YEAR, 2.37847173959480950e-03 * DAYS_PER_YEAR),
            vz: -2.96589568540237556e-05 * DAYS_PER_YEAR,
            mass: 4.36624404335156298e-05 * SOLAR_MASS,
        }
    }

    fn new_neptune() -> Self {
        Body {
            xy: f64x2::new(1.53796971148509165e+01, -2.59193146099879641e+01),
            z: 1.79258772950371181e-01,
            vxy: f64x2::new(2.68067772490389322e-03 * DAYS_PER_YEAR, 1.62824170038242295e-03 * DAYS_PER_YEAR),
            vz: -9.51592254519715870e-05 * DAYS_PER_YEAR,
            mass: 5.15138902046611451e-05 * SOLAR_MASS,
        }
    }

    fn new_sun_no_velocity() -> Self {
        Body {
            xy: f64x2::new(0., 0.),
            z: 0.,
            vxy: f64x2::new(0., 0.),
            vz: 0.,
            mass: SOLAR_MASS,
        }
    }
}
