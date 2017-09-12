extern crate simd;

use simd::x86::sse2::f64x2;
use std::env;
use std::f64::consts::PI;

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

struct Bodies(Body, Body, Body, Body, Body);

impl Bodies {
    fn new() -> Self {
        Bodies(Body::new_sun(), Body::new_jupiter(), Body::new_saturn(), Body::new_uranus(), Body::new_neptune())
    }

    fn get_energy(&self) -> f64 {
        let mut e = 0.;
        let bodies = [&self.0, &self.1, &self.2, &self.3, &self.4];
        for i in 0..bodies.len() {
            let body_i = bodies[i];
            let vxy_pow = body_i.vxy * body_i.vxy;
            e += 0.5 * body_i.mass * (
                vxy_pow.extract(0) +
                    vxy_pow.extract(1) +
                    body_i.vz.powi(2));
            for j in i + 1..bodies.len() {
                let body_j = bodies[j];
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
            update_v(&mut self.0, &mut self.1);
            update_v(&mut self.0, &mut self.2);
            update_v(&mut self.0, &mut self.3);
            update_v(&mut self.0, &mut self.4);
            update_v(&mut self.1, &mut self.2);
            update_v(&mut self.1, &mut self.3);
            update_v(&mut self.1, &mut self.4);
            update_v(&mut self.2, &mut self.3);
            update_v(&mut self.2, &mut self.4);
            update_v(&mut self.3, &mut self.4);
            update_pos(&mut self.0);
            update_pos(&mut self.1);
            update_pos(&mut self.2);
            update_pos(&mut self.3);
            update_pos(&mut self.4);
        }
    }
}

#[inline(always)]
fn update_v(body_i: &mut Body, body_j: &mut Body) {
    let dxy = body_i.xy - body_j.xy;
    let dz = body_i.z - body_j.z;

    let dxy_2 = dxy * dxy;
    let d = dxy_2.extract(0) + dxy_2.extract(1) + dz * dz;
    let mag = TIME_STEP / (d.sqrt() * d);

    let m_j_mag = body_j.mass * mag;
    body_i.vxy = body_i.vxy - dxy * f64x2::splat(m_j_mag);
    body_i.vz -= dz * m_j_mag;

    let m_i_mag = body_i.mass * mag;
    body_j.vxy = body_j.vxy + dxy * f64x2::splat(m_i_mag);
    body_j.vz += dz * m_i_mag;
}

#[inline(always)]
fn update_pos(body: &mut Body) {
    body.xy = body.xy + body.vxy * f64x2::splat(TIME_STEP);
    body.z += body.vz * TIME_STEP;
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

    fn new_sun() -> Self {
        Body {
            xy: f64x2::new(0., 0.),
            z: 0.,
            vxy: f64x2::new(-3.8766340719874267e-04, -3.2753590371765707e-03),
            vz: 2.3935734080003002e-05,
            mass: SOLAR_MASS,
        }
    }
}
