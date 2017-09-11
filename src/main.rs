use std::env;
use std::f64::consts::PI;

const BODIES: usize = 5;
const SOLAR_MASS: f64 = 4. * PI * PI;
const DAYS_PER_YEAR: f64 = 365.24;
const TIME_STEP: f64 = 0.01;

macro_rules! looper {
    ($i: expr, $($j: expr),+) => {{
        $(
            println!("{} : {}", $i, $j);
        )+
        looper!($($j),+)
    }};

    ($i: expr) => {{}}
}

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
            bodies[0].vx -= bodies[i].vx * bodies[i].mass / SOLAR_MASS;
            bodies[0].vy -= bodies[i].vy * bodies[i].mass / SOLAR_MASS;
            bodies[0].vz -= bodies[i].vz * bodies[i].mass / SOLAR_MASS;
        }
        Bodies { bodies }
    }

    fn get_energy(&self) -> f64 {
        let mut e = 0.;
        for i in 0..BODIES {
            let body_i = &self.bodies[i];
            e += 0.5 * body_i.mass * (
                body_i.vx.powi(2) +
                body_i.vy.powi(2) +
                body_i.vz.powi(2));
            for j in i + 1..BODIES {
                let body_j = &self.bodies[j];
                let d =
                    (body_i.x - body_j.x).powi(2) +
                    (body_i.y - body_j.y).powi(2) +
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
                    let dx = self.bodies[i].x - self.bodies[j].x;
                    let dy = self.bodies[i].y - self.bodies[j].y;
                    let dz = self.bodies[i].z - self.bodies[j].z;

                    let d = dx * dx + dy * dy + dz * dz;
                    let mag = TIME_STEP / (d.sqrt() * d);

                    self.bodies[i].vx -= dx * self.bodies[j].mass * mag;
                    self.bodies[i].vy -= dy * self.bodies[j].mass * mag;
                    self.bodies[i].vz -= dz * self.bodies[j].mass * mag;

                    self.bodies[j].vx += dx * self.bodies[i].mass * mag;
                    self.bodies[j].vy += dy * self.bodies[i].mass * mag;
                    self.bodies[j].vz += dz * self.bodies[i].mass * mag;
                }
                self.bodies[i].x += TIME_STEP * self.bodies[i].vx;
                self.bodies[i].y += TIME_STEP * self.bodies[i].vy;
                self.bodies[i].z += TIME_STEP * self.bodies[i].vz;
            }
        }

//        for _ in 0..steps {
//            for i in 0..BODIES {
//                let (bodies_i, bodies_j) = self.bodies.split_at_mut(i+1);
//                let body_i = &mut bodies_i[i];
//                for j in 0..BODIES - i - 1 {
//                    let body_j = &mut bodies_j[j];
//                    let dx = body_i.x - body_j.x;
//                    let dy = body_i.y - body_j.y;
//                    let dz = body_i.z - body_j.z;
//
//                    let d = dx * dx + dy * dy + dz * dz;
//                    let m = TIME_STEP / (d.sqrt() * d);
//
//                    body_i.vx -= dx * body_j.mass * m;
//                    body_i.vy -= dy * body_j.mass * m;
//                    body_i.vz -= dz * body_j.mass * m;
//
//                    body_j.vx += dx * body_i.mass * m;
//                    body_j.vy += dy * body_i.mass * m;
//                    body_j.vz += dz * body_i.mass * m;
//                }
//                body_i.x += TIME_STEP * body_i.vx;
//                body_i.y += TIME_STEP * body_i.vy;
//                body_i.z += TIME_STEP * body_i.vz;
//            }
//        }
    }
}

struct Body {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    mass: f64,
}

impl Body {
    fn new_jupiter() -> Self {
        Body {
            x: 4.84143144246472090e+00,
            y: -1.16032004402742839e+00,
            z: -1.03622044471123109e-01,
            vx: 1.66007664274403694e-03 * DAYS_PER_YEAR,
            vy: 7.69901118419740425e-03 * DAYS_PER_YEAR,
            vz: -6.90460016972063023e-05 * DAYS_PER_YEAR,
            mass: 9.54791938424326609e-04 * SOLAR_MASS,
        }
    }

    fn new_saturn() -> Self {
        Body {
            x: 8.34336671824457987e+00,
            y: 4.12479856412430479e+00,
            z: -4.03523417114321381e-01,
            vx: -2.76742510726862411e-03 * DAYS_PER_YEAR,
            vy: 4.99852801234917238e-03 * DAYS_PER_YEAR,
            vz: 2.30417297573763929e-05 * DAYS_PER_YEAR,
            mass: 2.85885980666130812e-04 * SOLAR_MASS,
        }
    }

    fn new_uranus() -> Self {
        Body {
            x: 1.28943695621391310e+01,
            y: -1.51111514016986312e+01,
            z: -2.23307578892655734e-01,
            vx: 2.96460137564761618e-03 * DAYS_PER_YEAR,
            vy: 2.37847173959480950e-03 * DAYS_PER_YEAR,
            vz: -2.96589568540237556e-05 * DAYS_PER_YEAR,
            mass: 4.36624404335156298e-05 * SOLAR_MASS,
        }
    }

    fn new_neptune() -> Self {
        Body {
            x: 1.53796971148509165e+01,
            y: -2.59193146099879641e+01,
            z: 1.79258772950371181e-01,
            vx: 2.68067772490389322e-03 * DAYS_PER_YEAR,
            vy: 1.62824170038242295e-03 * DAYS_PER_YEAR,
            vz: -9.51592254519715870e-05 * DAYS_PER_YEAR,
            mass: 5.15138902046611451e-05 * SOLAR_MASS,
        }
    }

    fn new_sun_no_velocity() -> Self {
        Body {
            x: 0.,
            y: 0.,
            z: 0.,
            vx: 0.,
            vy: 0.,
            vz: 0.,
            mass: SOLAR_MASS,
        }
    }
}
