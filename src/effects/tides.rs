use super::super::constants::{K2, PI};
use super::super::tools;
use super::super::Axes;
use super::super::Particle;
use super::EvolutionType;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleInputParameters {
    pub dissipation_factor: f64,
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub love_number: f64, // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the
    // susceptibility of its shape to change in response to a tidal potential.
    pub creep_coplanar_tides_input_parameters: CreepCoplanarTidesInputParticleParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct CreepCoplanarTidesInputParticleParameters {
    pub uniform_viscosity_coefficient: f64,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleInternalParameters {
    pub distance: f64,
    pub radial_velocity: f64,
    pub scaled_dissipation_factor: f64, // sigma (dissipation_factor_scale*dissipation_factor)
    pub scalar_product_of_vector_position_with_stellar_spin: f64,
    pub scalar_product_of_vector_position_with_planetary_spin: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_stellar_tide: f64,
    pub orthogonal_component_of_the_tidal_force_due_to_planetary_tide: f64,
    pub radial_component_of_the_tidal_force: f64,
    pub radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass: f64, // Needed to compute denergy_dt
    pub denergy_dt: f64, // Only for history output
    pub lag_angle: f64, // Used by EvolutionType::BolmontMathis2016, EvolutionType::GalletBolmont2017 and EvolutionType::LeconteChabrier2013(true)
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleParameters {
    pub input: TidesParticleInputParameters,
    pub internal: TidesParticleInternalParameters,
    pub output: TidesParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct TidesParticleCoordinates {
    // Positions/velocities in a heliocentric frame
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum TidesEffect {
    CentralBody,
    OrbitingBody,
    CreepCoplanarOrbitingBody,
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Tides {
    pub effect: TidesEffect,
    pub parameters: TidesParticleParameters,
    pub coordinates: TidesParticleCoordinates,
}

impl Tides {
    pub fn new(
        effect: TidesEffect,
        dissipation_factor: f64,
        dissipation_factor_scale: f64,
        love_number: f64,
        uniform_viscosity_coefficient: f64,
    ) -> Tides {
        Tides {
            effect: effect,
            parameters: TidesParticleParameters {
                input: TidesParticleInputParameters {
                    dissipation_factor: dissipation_factor,
                    dissipation_factor_scale: dissipation_factor_scale,
                    love_number: love_number,
                    // ////////////////////////////////////////////////////////////////////////////////
                    creep_coplanar_tides_input_parameters:
                        CreepCoplanarTidesInputParticleParameters {
                            uniform_viscosity_coefficient: uniform_viscosity_coefficient,
                        },
                    // ////////////////////////////////////////////////////////////////////////////////
                },
                internal: TidesParticleInternalParameters {
                    distance: 0.,
                    radial_velocity: 0.,
                    scaled_dissipation_factor: dissipation_factor_scale * dissipation_factor,
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    orthogonal_component_of_the_tidal_force_due_to_stellar_tide: 0.,
                    orthogonal_component_of_the_tidal_force_due_to_planetary_tide: 0.,
                    radial_component_of_the_tidal_force: 0.,
                    radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass:
                        0.,
                    denergy_dt: 0., // Only for history output
                    lag_angle: 0.,  // It will be initialized the first time the evolver is called
                },
                output: TidesParticleOutputParameters {
                    acceleration: Axes {
                        x: 0.,
                        y: 0.,
                        z: 0.,
                    },
                    dangular_momentum_dt: Axes {
                        x: 0.,
                        y: 0.,
                        z: 0.,
                    },
                },
            },
            coordinates: TidesParticleCoordinates {
                position: Axes {
                    x: 0.,
                    y: 0.,
                    z: 0.,
                },
                velocity: Axes {
                    x: 0.,
                    y: 0.,
                    z: 0.,
                },
            },
        }
    }
}

pub fn initialize(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let TidesEffect::CentralBody = host_particle.tides.effect {
        host_particle
            .tides
            .parameters
            .internal
            .scalar_product_of_vector_position_with_stellar_spin = 0.;
        host_particle
            .tides
            .parameters
            .internal
            .scalar_product_of_vector_position_with_planetary_spin = 0.;
        host_particle.tides.parameters.output.acceleration.x = 0.;
        host_particle.tides.parameters.output.acceleration.y = 0.;
        host_particle.tides.parameters.output.acceleration.z = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
        host_particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody = particle.tides.effect {
                //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {   // Check with the syntax of Alex for this
                //if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
                particle
                    .tides
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_stellar_spin =
                    particle.tides.coordinates.position.x * host_particle.spin.x
                        + particle.tides.coordinates.position.y * host_particle.spin.y
                        + particle.tides.coordinates.position.z * host_particle.spin.z;
                particle
                    .tides
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_planetary_spin =
                    particle.tides.coordinates.position.x * particle.spin.x
                        + particle.tides.coordinates.position.y * particle.spin.y
                        + particle.tides.coordinates.position.z * particle.spin.z;
                particle.tides.parameters.output.acceleration.x = 0.;
                particle.tides.parameters.output.acceleration.y = 0.;
                particle.tides.parameters.output.acceleration.z = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
            }
            if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
                //if let TidesEffect::ConstTimeLagOrbitingBody = particle.tides.effect {   // Check with the syntax of Alex for this
                //if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
                particle
                    .tides
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_stellar_spin =
                    particle.tides.coordinates.position.x * host_particle.spin.x
                        + particle.tides.coordinates.position.y * host_particle.spin.y
                        + particle.tides.coordinates.position.z * host_particle.spin.z;
                particle
                    .tides
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_planetary_spin =
                    particle.tides.coordinates.position.x * particle.spin.x
                        + particle.tides.coordinates.position.y * particle.spin.y
                        + particle.tides.coordinates.position.z * particle.spin.z;
                particle.tides.parameters.output.acceleration.x = 0.;
                particle.tides.parameters.output.acceleration.y = 0.;
                particle.tides.parameters.output.acceleration.z = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.x = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.y = 0.;
                particle.tides.parameters.output.dangular_momentum_dt.z = 0.;
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let TidesEffect::CentralBody = host_particle.tides.effect {
        // Inertial to Heliocentric positions/velocities
        host_particle.tides.coordinates.position.x = 0.;
        host_particle.tides.coordinates.position.y = 0.;
        host_particle.tides.coordinates.position.z = 0.;
        host_particle.tides.coordinates.velocity.x = 0.;
        host_particle.tides.coordinates.velocity.y = 0.;
        host_particle.tides.coordinates.velocity.z = 0.;
        host_particle.tides.parameters.internal.distance = 0.;
        host_particle.tides.parameters.internal.radial_velocity = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody = particle.tides.effect {
                particle.tides.coordinates.position.x =
                    particle.inertial_position.x - host_particle.inertial_position.x;
                particle.tides.coordinates.position.y =
                    particle.inertial_position.y - host_particle.inertial_position.y;
                particle.tides.coordinates.position.z =
                    particle.inertial_position.z - host_particle.inertial_position.z;
                particle.tides.coordinates.velocity.x =
                    particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.tides.coordinates.velocity.y =
                    particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.tides.coordinates.velocity.z =
                    particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.tides.parameters.internal.distance =
                    (particle.tides.coordinates.position.x.powi(2)
                        + particle.tides.coordinates.position.y.powi(2)
                        + particle.tides.coordinates.position.z.powi(2))
                    .sqrt();
                particle.tides.parameters.internal.radial_velocity =
                    (particle.tides.coordinates.position.x * particle.tides.coordinates.velocity.x
                        + particle.tides.coordinates.position.y
                            * particle.tides.coordinates.velocity.y
                        + particle.tides.coordinates.position.z
                            * particle.tides.coordinates.velocity.z)
                        / particle.tides.parameters.internal.distance;
            }
            if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
                particle.tides.coordinates.position.x =
                    particle.inertial_position.x - host_particle.inertial_position.x;
                particle.tides.coordinates.position.y =
                    particle.inertial_position.y - host_particle.inertial_position.y;
                particle.tides.coordinates.position.z =
                    particle.inertial_position.z - host_particle.inertial_position.z;
                particle.tides.coordinates.velocity.x =
                    particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.tides.coordinates.velocity.y =
                    particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.tides.coordinates.velocity.z =
                    particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.tides.parameters.internal.distance =
                    (particle.tides.coordinates.position.x.powi(2)
                        + particle.tides.coordinates.position.y.powi(2)
                        + particle.tides.coordinates.position.z.powi(2))
                    .sqrt();
                particle.tides.parameters.internal.radial_velocity =
                    (particle.tides.coordinates.position.x * particle.tides.coordinates.velocity.x
                        + particle.tides.coordinates.position.y
                            * particle.tides.coordinates.velocity.y
                        + particle.tides.coordinates.position.z
                            * particle.tides.coordinates.velocity.z)
                        / particle.tides.parameters.internal.distance;
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let TidesEffect::CentralBody = host_particle.tides.effect {
        host_particle.tides.coordinates.position = host_particle.heliocentric_position;
        host_particle.tides.coordinates.velocity = host_particle.heliocentric_velocity;
        host_particle.tides.parameters.internal.distance = host_particle.heliocentric_distance;
        host_particle.tides.parameters.internal.radial_velocity =
            host_particle.heliocentric_radial_velocity;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let TidesEffect::OrbitingBody = particle.tides.effect {
                particle.tides.coordinates.position = particle.heliocentric_position;
                particle.tides.coordinates.velocity = particle.heliocentric_velocity;
                particle.tides.parameters.internal.distance = particle.heliocentric_distance;
                particle.tides.parameters.internal.radial_velocity =
                    particle.heliocentric_radial_velocity;
            }
            if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
                particle.tides.coordinates.position = particle.heliocentric_position;
                particle.tides.coordinates.velocity = particle.heliocentric_velocity;
                particle.tides.parameters.internal.distance = particle.heliocentric_distance;
                particle.tides.parameters.internal.radial_velocity =
                    particle.heliocentric_radial_velocity;
            }
        }
    }
}

pub fn calculate_planet_dependent_dissipation_factors(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>,
) {
    match tidal_host_particle.evolution {
        EvolutionType::BolmontMathis2016(_)
        | EvolutionType::GalletBolmont2017(_)
        | EvolutionType::LeconteChabrier2013(true) => {
            let star_norm_spin_vector = tidal_host_particle.norm_spin_vector_2.sqrt();
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                //
                //// Excitation frequency needed by the model based on the
                // instantaneous frequency (using positions, velocities and spins)
                //let frequency = (particle.tides.coordinates.velocity.x - tidal_host_particle.spin.y*particle.tides.coordinates.position.z + tidal_host_particle.spin.z*particle.tides.coordinates.position.y).powi(2)
                //+ (particle.tides.coordinates.velocity.y - tidal_host_particle.spin.z*particle.tides.coordinates.position.x + tidal_host_particle.spin.x*particle.tides.coordinates.position.z).powi(2)
                //+ (particle.tides.coordinates.velocity.z - tidal_host_particle.spin.x*particle.tides.coordinates.position.y + tidal_host_particle.spin.y*particle.tides.coordinates.position.x).powi(2);
                //let inverse_of_half_the_excitation_frequency = particle.tides.parameters.internal.distance / frequency;
                // NOTE:  two_times_the_inverse_of_the_excitation_frequency: 2/w
                //        inverse_of_half_the_excitation_frequency : 1/(w/2)
                //
                //// Excitation frequency needed by the model based on the
                // mean frequency (using mean motion and spin).
                //
                // NOTE: The model is already here being used outside the
                // validity domain, it seems not justified to use an
                // instantaneous frequency.
                let gm = tidal_host_particle.mass_g + particle.mass_g;
                let (perihelion_distance, eccentricity) =
                    tools::calculate_perihelion_distance_and_eccentricity(
                        gm,
                        particle.tides.coordinates.position,
                        particle.tides.coordinates.velocity,
                    );
                let mean_motion =
                    gm.sqrt() * (perihelion_distance / (1.0 - eccentricity)).powf(-1.5);
                let half_the_excitation_frequency = (star_norm_spin_vector - mean_motion).abs();
                let inverse_of_half_the_excitation_frequency = 1. / half_the_excitation_frequency;

                let planet_dependent_dissipation_factor = tidal_host_particle
                    .tides
                    .parameters
                    .input
                    .dissipation_factor_scale
                    * 2.0
                    * K2
                    * tidal_host_particle.tides.parameters.internal.lag_angle
                    * inverse_of_half_the_excitation_frequency
                    / (3.0 * tidal_host_particle.radius.powi(5));

                star_planet_dependent_dissipation_factors
                    .insert(particle.id.clone(), planet_dependent_dissipation_factor);
                //println!("Insert {} in {}", planet_dependent_dissipation_factor, particle.id);
            }
            //panic!("Please, contact Posidonius authors before using BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true) evolutionary models. They may not be ready yet for scientific explotation.")
        }
        _ => {}
    }
}

pub fn planet_dependent_dissipation_factor(
    star_planet_dependent_dissipation_factors: &HashMap<usize, f64>,
    id: &usize,
    evolution: EvolutionType,
    scaled_dissipation_factor: f64,
) -> f64 {
    match evolution {
        EvolutionType::BolmontMathis2016(_)
        | EvolutionType::GalletBolmont2017(_)
        | EvolutionType::LeconteChabrier2013(true) => {
            match star_planet_dependent_dissipation_factors.get(id) {
                Some(&value) => value,
                _ => scaled_dissipation_factor, // This should not happen
            }
        }
        _ => scaled_dissipation_factor,
    }
}

#[allow(dead_code)]
pub fn calculate_true_anomaly(eccentricity: f64, mean_anomaly: f64) -> f64 {
    let mut difference: f64 = 10.0;
    let precision: f64 = 1.0e-13;
    let mut eccentric_anomaly_0: f64 = mean_anomaly;
    let mut eccentric_anomaly_1: f64;
    let true_anomaly: f64;

    while difference > precision {
        eccentric_anomaly_1 = mean_anomaly + eccentricity * eccentric_anomaly_0.sin();
        difference = (eccentric_anomaly_1 - eccentric_anomaly_0).abs();
        eccentric_anomaly_0 = eccentric_anomaly_1;
    }

    //   println!("{}", mean_anomaly);
    //   println!("{}", eccentric_anomaly_0);

    true_anomaly = 2.0
        * (((1.0 + eccentricity) / (1.0 - eccentricity)).powf(1.0 / 2.0)
            * (eccentric_anomaly_0 / 2.0).tan())
        .atan();

    //   println!("{}", true_anomaly);

    return true_anomaly;
}

#[allow(dead_code)]
pub fn calculate_cayley_coefficients(e: f64) -> [[f64; 15]; 2] {
    let mut cayley = [[0.; 15]; 2];

    let e2 = e.powi(2);
    let e3 = e * e2;
    let e4 = e2 * e2;
    let e5 = e4 * e;
    let e6 = e5 * e;
    let e7 = e6 * e;

    // Coefficients E_(0,k) of FM2015 arvix

    cayley[0][0] = 432091.0 / 30720.0 * e7;
    cayley[0][1] = 3167.0 / 320.0 * e6;
    cayley[0][2] = 1773.0 / 256.0 * e5 - 4987.0 / 6144.0 * e7;
    cayley[0][3] = 77.0 / 16.0 * e4 + 129.0 / 160.0 * e6;
    cayley[0][4] = 53.0 / 16.0 * e3 + 393.0 / 256.0 * e5 + 24753.0 / 10240.0 * e7;
    cayley[0][5] = 9.0 / 4.0 * e2 + 7.0 / 4.0 * e4 + 141.0 / 64.0 * e6;
    cayley[0][6] = 3.0 / 2.0 * e + 27.0 / 16.0 * e3 - 261.0 / 128.0 * e5 + 14309.0 / 6144.0 * e7;
    cayley[0][7] = 1.0 + 3.0 / 2.0 * e2 + 15.0 / 8.0 * e4 - 35.0 / 16.0 * e6;
    cayley[0][8] = cayley[0][6];
    cayley[0][9] = cayley[0][5];
    cayley[0][10] = cayley[0][4];
    cayley[0][11] = cayley[0][3];
    cayley[0][12] = cayley[0][2];
    cayley[0][13] = cayley[0][1];
    cayley[0][14] = cayley[0][0];

    // Coefficients E_(2,k) of FM2015 arvix

    cayley[1][0] = 12144273.0 / 71680.0 * e7;
    cayley[1][1] = 73369.0 / 720.0 * e6;
    cayley[1][2] = 228347.0 / 3840.0 * e5 - 3071075.0 / 18432.0 * e7;
    cayley[1][3] = 533.0 / 16.0 * e4 - 13827.0 / 160.0 * e6;
    cayley[1][4] = 845.0 / 48.0 * e3 - 32525.0 / 768.0 * e5 + 208225.0 / 6144.0 * e7;
    cayley[1][5] = 17.0 / 2.0 * e2 - 115.0 / 6.0 * e4 + 601.0 / 48.0 * e6;
    cayley[1][6] = 7.0 / 2.0 * e - 123.0 / 16.0 * e3 + 489.0 / 128.0 * e5 - 1763.0 / 2048.0 * e7;
    cayley[1][7] = 1.0 - 5.0 / 2.0 * e2 + 13.0 / 16.0 * e4 - 35.0 / 288.0 * e6;
    cayley[1][8] = -1.0 / 2.0 * e + 1.0 / 16.0 * e3 - 5.0 / 384.0 * e5 - 143.0 / 18432.0 * e7;
    cayley[1][9] = 0.0;
    cayley[1][10] = 1.0 / 48.0 * e3 + 11.0 / 768.0 * e5 + 313.0 / 30720.0 * e7;
    cayley[1][11] = 1.0 / 24.0 * e4 + 7.0 / 240.0 * e6;
    cayley[1][12] = 81.0 / 1280.0 * e5 + 81.0 / 2048.0 * e7;
    cayley[1][13] = 4.0 / 45.0 * e6;
    cayley[1][14] = 15625.0 / 129024.0 * e7;

    return cayley;
}

#[allow(dead_code)]
pub fn calculate_particle_shape(
    semimajor_axis: f64,
    eccentricity: f64,
    mean_anomaly: f64,
    orbital_period: f64,
    primary_mass: f64,
    companion_mass: f64,
    radius: f64,
    uniform_viscosity_coefficient: f64,
    primary_spin: f64,
) -> Axes {
    let epsilon_rho_bar: f64;
    let epsilon_z_bar: f64;
    let mean_motion: f64;
    let relaxation_factor: f64;
    let (
        x_component_of_the_equatorial_prolateness,
        y_component_of_the_equatorial_prolateness,
        epsilon_z,
    ): (f64, f64, f64);
    let alpha: f64;
    let mean_motion_to_relaxation_factor_ratio: f64;
    let (
        equatorial_prolateness_constant_coefficient,
        equatorial_prolateness_cosine_coefficient,
        equatorial_prolateness_sine_coefficient,
    ): (f64, f64, f64);
    let (
        polar_oblateness_constant_coefficient,
        polar_oblateness_cosine_coefficient,
        polar_oblateness_sine_coefficient,
    ): (f64, f64, f64);
    let (delta_constant_coefficient, delta_cosine_coefficient, delta_sine_coefficient): (
        f64,
        f64,
        f64,
    );
    let (ep_rho, delta, delta2, epsilon_z): (f64, f64, f64, f64);
    let companion_mass_to_sum_of_masses_ratio: f64;
    let pseudo_synchronization_frequency: f64;

    relaxation_factor = 3.0 * K2 * primary_mass * primary_mass
        / (8.0 * PI * radius.powi(4) * uniform_viscosity_coefficient);
    mean_motion = 2.0 * PI / orbital_period;
    epsilon_rho_bar =
        15.0 * companion_mass * radius.powi(3) / (4.0 * primary_mass * semimajor_axis.powi(3));
    epsilon_z_bar = 5.0 * primary_spin * primary_spin * radius.powi(3) / (4.0 * K2 * primary_mass);
    companion_mass_to_sum_of_masses_ratio = companion_mass / (primary_mass + companion_mass);

         // This part is important, adjust properly!

    if relaxation_factor < mean_motion {
        pseudo_synchronization_frequency = 0.9 * mean_motion; // 1.0 * mean_motion + 6.0 * eccentricity * eccentricity * mean_motion;
    } else {
        pseudo_synchronization_frequency =
            1.0 * mean_motion + 7.0 * eccentricity * eccentricity * mean_motion; // 1.0 * mean_motion + 8.0 * eccentricity * eccentricity * mean_motion;
    }

    if primary_spin < pseudo_synchronization_frequency {
        // be careful here
        // The Epsilons are taken from Equation 49 of Folonier et al. (2018)

        alpha = 1.0 - (3.0 * companion_mass_to_sum_of_masses_ratio * epsilon_rho_bar);
        mean_motion_to_relaxation_factor_ratio = mean_motion / relaxation_factor;
        // mean_motion_to_relaxation_factor_ratio_2 = mean_motion_to_relaxation_factor_ratio.powi(2);
        // eccentricity_2 = eccentricity.powi(2);
        // alpha_2 = alpha.powi(2);

        equatorial_prolateness_constant_coefficient = epsilon_rho_bar
            * (1.0 + 1.5 * eccentricity.powi(2)
                - (4.0 * mean_motion_to_relaxation_factor_ratio.powi(2) * eccentricity.powi(2)
                    / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2))));
        equatorial_prolateness_cosine_coefficient = 3.0 * epsilon_rho_bar * eccentricity
            / (1.0 + mean_motion_to_relaxation_factor_ratio.powi(2));
        equatorial_prolateness_sine_coefficient =
            3.0 * epsilon_rho_bar * eccentricity * mean_motion_to_relaxation_factor_ratio
                / (1.0 + mean_motion_to_relaxation_factor_ratio.powi(2));

        polar_oblateness_constant_coefficient = epsilon_rho_bar
            * (0.5 + 0.75 * eccentricity.powi(2))
            + epsilon_z_bar
                * (1.0
                    + 12.0
                        * eccentricity.powi(2)
                        * (1.0 + alpha * mean_motion_to_relaxation_factor_ratio.powi(2))
                        / ((1.0 + mean_motion_to_relaxation_factor_ratio.powi(2))
                            * (1.0
                                + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2)))
                    + 2.0
                        * (1.0 - alpha)
                        * (1.0 - alpha)
                        * mean_motion_to_relaxation_factor_ratio.powi(2)
                        * eccentricity.powi(2)
                        / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2)));
        polar_oblateness_cosine_coefficient = (1.5 * epsilon_rho_bar * eccentricity
            / (1.0 + mean_motion_to_relaxation_factor_ratio.powi(2)))
            * (1.0
                - (16.0
                    * companion_mass_to_sum_of_masses_ratio
                    * mean_motion_to_relaxation_factor_ratio.powi(2)
                    * epsilon_z_bar
                    / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2))));
        polar_oblateness_sine_coefficient = (1.5 * epsilon_rho_bar * eccentricity
            / (1.0 + mean_motion_to_relaxation_factor_ratio.powi(2)))
            * (mean_motion_to_relaxation_factor_ratio
                + (8.0
                    * companion_mass_to_sum_of_masses_ratio
                    * mean_motion_to_relaxation_factor_ratio
                    * (1.0 - alpha * mean_motion_to_relaxation_factor_ratio.powi(2))
                    * epsilon_z_bar
                    / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2))));

        delta_constant_coefficient =
            (3.0 * mean_motion_to_relaxation_factor_ratio * eccentricity.powi(2)
                / (1.0 + mean_motion_to_relaxation_factor_ratio.powi(2)))
                * (2.0 + (1.0 + alpha) * mean_motion_to_relaxation_factor_ratio.powi(2))
                / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2));
        delta_cosine_coefficient = -2.0 * mean_motion_to_relaxation_factor_ratio * eccentricity
            / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2));
        delta_sine_coefficient =
            -2.0 * eccentricity * mean_motion_to_relaxation_factor_ratio.powi(2) * alpha
                / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio.powi(2));

        ep_rho = equatorial_prolateness_constant_coefficient
            + equatorial_prolateness_cosine_coefficient * mean_anomaly.cos()
            + equatorial_prolateness_sine_coefficient * mean_anomaly.sin();

        epsilon_z = polar_oblateness_constant_coefficient
            + polar_oblateness_cosine_coefficient * mean_anomaly.cos()
            + polar_oblateness_sine_coefficient * mean_anomaly.sin();

        delta = delta_constant_coefficient
            + delta_cosine_coefficient * mean_anomaly.cos()
            + delta_sine_coefficient * mean_anomaly.sin();

        delta2 = 2.0 * delta;
        x_component_of_the_equatorial_prolateness = ep_rho * delta2.cos();
        y_component_of_the_equatorial_prolateness = ep_rho * delta2.sin();

        Axes {
            x: x_component_of_the_equatorial_prolateness,
            y: y_component_of_the_equatorial_prolateness,
            z: epsilon_z,
        }
    } else {
        let cayley = calculate_cayley_coefficients(eccentricity);
        let mut true_anomaly = calculate_true_anomaly(eccentricity, mean_anomaly);
        if true_anomaly < 0.0 {
            true_anomaly = true_anomaly + 2.0 * PI;
        }

        let mut sum_x_component_of_the_equatorial_prolateness: f64;
        let mut sum_y_component_of_the_equatorial_prolateness: f64;
        let mut sum_epsilon_z: f64;
        let mut equatorial_angles: f64;
        let mut polar_angles: f64;

        sum_x_component_of_the_equatorial_prolateness = 0.0;
        sum_y_component_of_the_equatorial_prolateness = 0.0;
        sum_epsilon_z = 0.0;

        let mut k: usize = 0;

        // for (k,factor) in cayley[1].iter().enumerate()

        while k < 15 {
            // remove

            let dk = k as f64;

            equatorial_angles = (9.0 - dk) * mean_anomaly - 2.0 * true_anomaly;
            polar_angles = (-7.0 + dk) * mean_anomaly;

            // common block (relaxation_factor * cayley[1][k] / (relaxation_factor.powi(2) + (2.0*primary_spin-(9.0-dk)*mean_motion).powi(2)))

             sum_x_component_of_the_equatorial_prolateness += relaxation_factor * cayley[1][k]
                / (relaxation_factor.powi(2)
                    + (2.0 * primary_spin - (9.0 - dk) * mean_motion).powi(2))
                * (relaxation_factor * (equatorial_angles).cos()
                    - (2.0 * primary_spin - (9.0 - dk) * mean_motion) * (equatorial_angles).sin());
             sum_y_component_of_the_equatorial_prolateness += relaxation_factor * cayley[1][k]
                / (relaxation_factor.powi(2)
                    + (2.0 * primary_spin - (9.0 - dk) * mean_motion).powi(2))
                * (relaxation_factor * (equatorial_angles).sin()
                    + (2.0 * primary_spin - (9.0 - dk) * mean_motion) * (equatorial_angles).cos());
             sum_epsilon_z += cayley[0][k]
                / ((dk - 7.0).powi(2) * mean_motion.powi(2) + relaxation_factor.powi(2))
                * ((dk - 7.0).powi(2) * mean_motion.powi(2) * (polar_angles).sin()
                    + relaxation_factor.powi(2) * (polar_angles).cos());
            //sum_x_component_of_the_equatorial_prolateness += 0.0;
            //sum_y_component_of_the_equatorial_prolateness += 0.0;
            //sum_epsilon_z += 0.0;
             

            // println!("{} {} {}", mean_anomaly * 180.0 / PI , relaxation_factor * cayley[1][k as usize] / (relaxation_factor.powi(2) + (2.0*primary_spin-(9.0-dk)*mean_motion).powi(2)) * (relaxation_factor * (equatorial_angles).cos() - (2.0*primary_spin-(9.0-dk)*mean_motion) * (equatorial_angles).sin()) , relaxation_factor * cayley[1][k as usize] / (relaxation_factor.powi(2) + (2.0*primary_spin-(9.0-dk)*mean_motion).powi(2)) * (relaxation_factor * (equatorial_angles).sin() + (2.0*primary_spin-(9.0-dk)*mean_motion) * (equatorial_angles).cos()));

            k += 1;
        }

        // println!("{}", 0.5 * sum_y_component_of_the_equatorial_prolateness / sum_x_component_of_the_equatorial_prolateness * 180.0 / PI);

        sum_epsilon_z = (sum_epsilon_z * epsilon_rho_bar / 2.0) + epsilon_z_bar; // In case tides + rotation are taken into account

        //sum_epsilon_z = epsilon_z_bar; // In case only rotational flattenings are taken into account

        x_component_of_the_equatorial_prolateness =
            sum_x_component_of_the_equatorial_prolateness * epsilon_rho_bar;
        y_component_of_the_equatorial_prolateness =
            sum_y_component_of_the_equatorial_prolateness * epsilon_rho_bar;
        epsilon_z = sum_epsilon_z;

        Axes {
            x: x_component_of_the_equatorial_prolateness,
            y: y_component_of_the_equatorial_prolateness,
            z: epsilon_z,
        }
    }
}

//#[allow(dead_code)]
//fn open_file() {
//    let path = Path::new("write_TTV_k2265b.txt");
//    let display = path.display();

// Open a file in write-only mode, returns `io::Result<File>`
//    let mut file = match File::create(&path) {
//        Err(why) => panic!("couldn't create {}: {}", display, why),
//        Ok(file) => file,
//    };
//}

//////////////////////////////////////////////////////////////////////////////
//// TIDES
pub fn calculate_torque_due_to_tides(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    central_body: bool,
) {
    let mut dangular_momentum_dt = Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let mut reference_spin = tidal_host_particle.spin.clone();
    let mut orthogonal_component_of_the_tidal_force: f64;
    let mut reference_rscalspin: f64;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            //if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
            if !central_body {
                reference_spin = particle.spin.clone();
                reference_rscalspin = particle
                    .tides
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_planetary_spin;
                orthogonal_component_of_the_tidal_force = particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
            } else {
                reference_rscalspin = particle
                    .tides
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_stellar_spin;
                orthogonal_component_of_the_tidal_force = particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
            }

            // distance to star
            let distance = particle.tides.parameters.internal.distance;

            //// Torque calculation (star)
            // - Equation 8-9 from Bolmont et al. 2015
            let torque_due_to_tides_x: f64 = orthogonal_component_of_the_tidal_force
                * (distance * reference_spin.x
                    - reference_rscalspin * particle.tides.coordinates.position.x / distance
                    - 1.0 / distance
                        * (particle.tides.coordinates.position.y
                            * particle.tides.coordinates.velocity.z
                            - particle.tides.coordinates.position.z
                                * particle.tides.coordinates.velocity.y));

            let torque_due_to_tides_y: f64 = orthogonal_component_of_the_tidal_force
                * (distance * reference_spin.y
                    - reference_rscalspin * particle.tides.coordinates.position.y / distance
                    - 1.0 / distance
                        * (particle.tides.coordinates.position.z
                            * particle.tides.coordinates.velocity.x
                            - particle.tides.coordinates.position.x
                                * particle.tides.coordinates.velocity.z));

            let torque_due_to_tidez_z: f64 = orthogonal_component_of_the_tidal_force
                * (distance * reference_spin.z
                    - reference_rscalspin * particle.tides.coordinates.position.z / distance
                    - 1.0 / distance
                        * (particle.tides.coordinates.position.x
                            * particle.tides.coordinates.velocity.y
                            - particle.tides.coordinates.position.y
                                * particle.tides.coordinates.velocity.x));

            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += factor * torque_due_to_tidez_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x =
                    factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y =
                    factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z =
                    factor * torque_due_to_tidez_z;
            }
        }
        if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
            let gm = tidal_host_particle.mass_g + particle.mass_g;
            let (
                semimajor_axis,
                _perihelion_distance,
                eccentricity,
                _inclination,
                _perihelion_longitude,
                _longitude_of_ascending_node,
                mean_anomaly,
                orbital_period,
            ) = tools::calculate_keplerian_orbital_elements(
                gm,
                particle.tides.coordinates.position,
                particle.tides.coordinates.velocity,
            );
            let spin = particle.spin.z;
            let distance = particle.tides.parameters.internal.distance;
            let true_anomaly = calculate_true_anomaly(eccentricity, mean_anomaly);

            //     println!("{} {}", perihelion_longitude, true_anomaly);

            let particle_shape = calculate_particle_shape(
                semimajor_axis,
                eccentricity,
                mean_anomaly,
                orbital_period,
                particle.mass,
                tidal_host_particle.mass,
                particle.radius,
                particle
                    .tides
                    .parameters
                    .input
                    .creep_coplanar_tides_input_parameters
                    .uniform_viscosity_coefficient,
                spin,
            );

            let torque_due_to_tides_x: f64 = 0.;
            let torque_due_to_tides_y: f64 = 0.;
            //let torque_due_to_tides_z: f64 = 0.0;

            // Torque expression for studying creep tide tidal despinning (use torque = 0 as above to study stat. rotation, makes code run faster)

            let torque_due_to_tides_z: f64 = 3.0 / 5.0
                * K2
                * particle.mass
                * tidal_host_particle.mass
                * particle.radius
                * particle.radius
                * particle_shape.y
                / distance.powi(3);
            
            //println!("{}", spin / (2.0 * PI / orbital_period));

            let factor = -1.0;
            if central_body {
                // Integration of the spin (total torque tides):
                dangular_momentum_dt.x += 0.;//factor * torque_due_to_tides_x;
                dangular_momentum_dt.y += 0.;//factor * torque_due_to_tides_y;
                dangular_momentum_dt.z += 0.;//factor * torque_due_to_tides_z;
            } else {
                particle.tides.parameters.output.dangular_momentum_dt.x =
                    factor * torque_due_to_tides_x;
                particle.tides.parameters.output.dangular_momentum_dt.y =
                    factor * torque_due_to_tides_y;
                particle.tides.parameters.output.dangular_momentum_dt.z =
                    factor * torque_due_to_tides_z;
            }
        }
    }

    if central_body {
        // - Equation 25 from Bolmont et al. 2015
        tidal_host_particle
            .tides
            .parameters
            .output
            .dangular_momentum_dt
            .x = dangular_momentum_dt.x;
        tidal_host_particle
            .tides
            .parameters
            .output
            .dangular_momentum_dt
            .y = dangular_momentum_dt.y;
        tidal_host_particle
            .tides
            .parameters
            .output
            .dangular_momentum_dt
            .z = dangular_momentum_dt.z;
    }
}

pub fn calculate_orthogonal_component_of_the_tidal_force(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>,
) {
    let mut tidal_host_particle = tidal_host_particle;
    let mut particles = particles;
    let mut more_particles = more_particles;
    let mut star_planet_dependent_dissipation_factors = star_planet_dependent_dissipation_factors;
    let central_body = true;
    calculate_orthogonal_component_of_the_tidal_force_for(
        central_body,
        &mut tidal_host_particle,
        &mut particles,
        &mut more_particles,
        &mut star_planet_dependent_dissipation_factors,
    );
    calculate_orthogonal_component_of_the_tidal_force_for(
        !central_body,
        &mut tidal_host_particle,
        &mut particles,
        &mut more_particles,
        &mut star_planet_dependent_dissipation_factors,
    );
}

fn calculate_orthogonal_component_of_the_tidal_force_for(
    central_body: bool,
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>,
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            //// Only calculate tides if planet is not in disk
            //if particle.disk_interaction_time == 0.0 {

            // (distance to star)^7
            let distance_7 = particle.tides.parameters.internal.distance.powi(7);

            //// Tidal force calculation (star) :: Only orthogonal component is needed
            if central_body {
                // - Third line of Equation 5 from Bolmont et al. 2015
                //   This expression has R**10 (instead of R**5 in Eq. 5)
                //   because it uses sigma (i.e., scaled_dissipation_factor)
                //   and not k2$\Delta$t (between k2$\Delta$t and sigma
                //   there is a R**5 factor as shown in Equation 28)
                //   - k2 is love number
                let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(
                    &star_planet_dependent_dissipation_factors,
                    &tidal_host_particle.id,
                    tidal_host_particle.evolution,
                    tidal_host_particle
                        .tides
                        .parameters
                        .internal
                        .scaled_dissipation_factor,
                );
                particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5
                    * (particle.mass.powi(2))
                    * (tidal_host_particle.radius.powi(10))
                    * star_scaled_dissipation_factor
                    / distance_7;
            } else {
                // - Second line of Equation 5 from Bolmont et al. 2015
                //   This expression has R**10 (instead of R**5 in Eq. 5)
                //   because it uses sigma (i.e., scaled_dissipation_factor)
                //   and not k2$\Delta$t (between k2$\Delta$t and sigma
                //   there is a R**5 factor as shown in Equation 28)
                //   - k2 is love number
                particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5
                    * (tidal_host_particle.mass.powi(2))
                    * (particle.radius.powi(10))
                    * particle.tides.parameters.internal.scaled_dissipation_factor
                    / distance_7;

                // SBC
                //println!("> {:e} {:e} {:e} {:e}", tidal_host_particle.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
                //println!("> {:e} {:e} {:e}", particle.tides.coordinates.position.x, particle.tides.coordinates.position.y, particle.tides.coordinates.position.z);
            }
            //} else {
            //if central_body {
            //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 0.0
            //} else{
            //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0
            //}
            //}
        }
    }
}

pub fn calculate_radial_component_of_the_tidal_force(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    star_planet_dependent_dissipation_factors: &mut HashMap<usize, f64>,
) {
    let star_mass_2 = tidal_host_particle.mass * tidal_host_particle.mass;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            let planet_mass_2 = particle.mass * particle.mass;
            // Conservative part of the radial tidal force
            let radial_component_of_the_tidal_force_conservative_part = -3.0 * K2
                / particle.tides.parameters.internal.distance.powi(7)
                * (planet_mass_2
                    * tidal_host_particle.radius.powi(5)
                    * tidal_host_particle.tides.parameters.input.love_number
                    + star_mass_2
                        * particle.radius.powi(5)
                        * particle.tides.parameters.input.love_number);

            // Dissipative part of the radial tidal force:
            let factor1 = -13.5 * particle.tides.parameters.internal.radial_velocity
                / particle.tides.parameters.internal.distance.powi(8);
            let star_scaled_dissipation_factor = planet_dependent_dissipation_factor(
                &star_planet_dependent_dissipation_factors,
                &particle.id,
                tidal_host_particle.evolution,
                tidal_host_particle
                    .tides
                    .parameters
                    .internal
                    .scaled_dissipation_factor,
            );
            let term1 = planet_mass_2
                * tidal_host_particle.radius.powi(10)
                * star_scaled_dissipation_factor;
            let term2 = star_mass_2
                * particle.radius.powi(10)
                * particle.tides.parameters.internal.scaled_dissipation_factor;
            // If we consider the star as a point mass (used for denergy_dt calculation):
            particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass =
                factor1 * term2;
            let radial_component_of_the_tidal_force_dissipative_part = particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass
                + factor1 * term1;

            // Sum of the dissipative and conservative part of the radial force
            // - First line Equation 5 from Bolmont et al. 2015
            particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force =
                radial_component_of_the_tidal_force_conservative_part
                    + radial_component_of_the_tidal_force_dissipative_part;
        }
    }
}

pub fn calculate_denergy_dt(particles: &mut [Particle], more_particles: &mut [Particle]) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            // - Equation 32 from Bolmont et al. 2015
            //// Instantaneous energy loss dE/dt due to tides
            //// in Msun.AU^2.day^(-3)
            //radial_tidal_force_for_energy_loss_calculation = factor1 * term2; // Ftidr_diss
            let factor2 = particle
                .tides
                .parameters
                .internal
                .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
                / particle.tides.parameters.internal.distance;
            particle.tides.parameters.internal.denergy_dt = -((1.0
                / particle.tides.parameters.internal.distance
                * (particle
                    .tides
                    .parameters
                    .internal
                    .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass
                    + factor2 * particle.tides.parameters.internal.radial_velocity))
                * (particle.tides.coordinates.position.x * particle.tides.coordinates.velocity.x
                    + particle.tides.coordinates.position.y
                        * particle.tides.coordinates.velocity.y
                    + particle.tides.coordinates.position.z
                        * particle.tides.coordinates.velocity.z)
                + factor2
                    * ((particle.spin.y * particle.tides.coordinates.position.z
                        - particle.spin.z * particle.tides.coordinates.position.y
                        - particle.tides.coordinates.velocity.x)
                        * particle.tides.coordinates.velocity.x
                        + (particle.spin.z * particle.tides.coordinates.position.x
                            - particle.spin.x * particle.tides.coordinates.position.z
                            - particle.tides.coordinates.velocity.y)
                            * particle.tides.coordinates.velocity.y
                        + (particle.spin.x * particle.tides.coordinates.position.y
                            - particle.spin.y * particle.tides.coordinates.position.x
                            - particle.tides.coordinates.velocity.z)
                            * particle.tides.coordinates.velocity.z))
                - (particle.tides.parameters.output.dangular_momentum_dt.x * particle.spin.x
                    + particle.tides.parameters.output.dangular_momentum_dt.y * particle.spin.y
                    + particle.tides.parameters.output.dangular_momentum_dt.z * particle.spin.z);
        }
    }
}

pub fn calculate_tidal_acceleration(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    let factor2 = 1. / tidal_host_particle.mass;
    let mut sum_total_tidal_force = Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody = particle.tides.effect {
            let factor1 = 1. / particle.mass;

            // - Equation 6 from Bolmont et al. 2015
            let factor3 = particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force
                + (particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
                    + particle
                        .tides
                        .parameters
                        .internal
                        .orthogonal_component_of_the_tidal_force_due_to_planetary_tide)
                    * particle.tides.parameters.internal.radial_velocity
                    / particle.tides.parameters.internal.distance;
            let total_tidal_force_x = factor3 * particle.tides.coordinates.position.x
                / particle.tides.parameters.internal.distance
                + particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
                    / particle.tides.parameters.internal.distance
                    * (tidal_host_particle.spin.y * particle.tides.coordinates.position.z
                        - tidal_host_particle.spin.z * particle.tides.coordinates.position.y
                        - particle.tides.coordinates.velocity.x)
                + particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
                    / particle.tides.parameters.internal.distance
                    * (particle.spin.y * particle.tides.coordinates.position.z
                        - particle.spin.z * particle.tides.coordinates.position.y
                        - particle.tides.coordinates.velocity.x);
            let total_tidal_force_y = factor3 * particle.tides.coordinates.position.y
                / particle.tides.parameters.internal.distance
                + particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
                    / particle.tides.parameters.internal.distance
                    * (tidal_host_particle.spin.z * particle.tides.coordinates.position.x
                        - tidal_host_particle.spin.x * particle.tides.coordinates.position.z
                        - particle.tides.coordinates.velocity.y)
                + particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
                    / particle.tides.parameters.internal.distance
                    * (particle.spin.z * particle.tides.coordinates.position.x
                        - particle.spin.x * particle.tides.coordinates.position.z
                        - particle.tides.coordinates.velocity.y);
            let total_tidal_force_z = factor3 * particle.tides.coordinates.position.z
                / particle.tides.parameters.internal.distance
                + particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
                    / particle.tides.parameters.internal.distance
                    * (tidal_host_particle.spin.x * particle.tides.coordinates.position.y
                        - tidal_host_particle.spin.y * particle.tides.coordinates.position.x
                        - particle.tides.coordinates.velocity.z)
                + particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
                    / particle.tides.parameters.internal.distance
                    * (particle.spin.x * particle.tides.coordinates.position.y
                        - particle.spin.y * particle.tides.coordinates.position.x
                        - particle.tides.coordinates.velocity.z);
            //println!("factor3 {:e} {:e} {:e}", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
            //println!("d {:e} vrad {:e}", particle.tides.parameters.internal., particle.tides.parameters.internal.radial_velocity);
            //println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);

            sum_total_tidal_force.x -= total_tidal_force_x; // CHANGE SIGN
            sum_total_tidal_force.y -= total_tidal_force_y; // CHANGE SIGN
            sum_total_tidal_force.z -= total_tidal_force_z; // CHANGE SIGN

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x;
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
            // println!("don't go here");
        }
        if let TidesEffect::CreepCoplanarOrbitingBody = particle.tides.effect {
            let factor1 = 1. / particle.mass;
            let gm = tidal_host_particle.mass_g + particle.mass_g;
            let (
                semimajor_axis,
                _perihelion_distance,
                eccentricity,
                _inclination,
                _perihelion_longitude,
                _longitude_of_ascending_node,
                mean_anomaly,
                orbital_period,
            ) = tools::calculate_keplerian_orbital_elements(
                gm,
                particle.tides.coordinates.position,
                particle.tides.coordinates.velocity,
            );
            let spin = particle.spin.z;
            let distance = particle.tides.parameters.internal.distance;
            let distance_4 = distance.powi(4);

            let particle_shape = calculate_particle_shape(
                semimajor_axis,
                eccentricity,
                mean_anomaly,
                orbital_period,
                particle.mass,
                tidal_host_particle.mass,
                particle.radius,
                particle
                    .tides
                    .parameters
                    .input
                    .creep_coplanar_tides_input_parameters
                    .uniform_viscosity_coefficient,
                spin,
            );

            // Factorize common terms;

           // Expression for tides + rotational flattenings

            let total_tidal_force_x = particle.tides.coordinates.position.x / distance
                * (-0.9
                    * K2
                    * particle.mass
                    * tidal_host_particle.mass
                    * particle.radius
                    * particle.radius
                    / distance_4
                    * particle_shape.x
                    - (3.0
                        * K2
                        * particle.mass
                        * tidal_host_particle.mass
                        * particle.radius
                        * particle.radius
                        / (5.0 * distance_4)
                        * particle_shape.z))
                - (3.0
                    * K2
                    * particle.tides.coordinates.position.y
                    * particle.mass
                    * tidal_host_particle.mass
                    * particle.radius
                    * particle.radius
                    / (5.0 * distance_4 * distance)
                    * particle_shape.y);

           // Expression for ONLY rotational flattenings

           // let total_tidal_force_x = particle.tides.coordinates.position.x / distance
           //         * (-3.0
           //             * K2
           //             * particle.mass
           //             * tidal_host_particle.mass
           //             * particle.radius
           //             * particle.radius
           //             / (5.0 * distance_4)
           //             * particle_shape.z);

           // Expression for tides + rotational flattenings                        

            let total_tidal_force_y = particle.tides.coordinates.position.y / distance
                * (-0.9
                    * K2
                    * particle.mass
                    * tidal_host_particle.mass
                    * particle.radius
                    * particle.radius
                    / distance_4
                    * particle_shape.x
                    - (3.0
                        * K2
                        * particle.mass
                        * tidal_host_particle.mass
                        * particle.radius
                        * particle.radius
                        / (5.0 * distance_4)
                        * particle_shape.z))
                + (3.0
                    * K2
                    * particle.tides.coordinates.position.x
                    * particle.mass
                    * tidal_host_particle.mass
                    * particle.radius
                    * particle.radius
                    / (5.0 * distance_4 * distance)
                    * particle_shape.y);

           // Expression for ONLY rotational flattenings                        

           // let total_tidal_force_y = particle.tides.coordinates.position.y / distance
           //     * (-3.0
           //         * K2
           //         * particle.mass
           //         * tidal_host_particle.mass
           //         * particle.radius
           //         * particle.radius
           //         / (5.0 * distance_4)
           //         * particle_shape.z);
           
            let total_tidal_force_z = 0.0;
            // println!("{} {} {}", particle_shape.x, particle_shape.y, particle_shape.z);
            //println!("{}", total_tidal_force_y);

            sum_total_tidal_force.x += total_tidal_force_x;
            sum_total_tidal_force.y += total_tidal_force_y;
            sum_total_tidal_force.z += total_tidal_force_z;

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.tides.parameters.output.acceleration.x = factor1 * total_tidal_force_x;
            particle.tides.parameters.output.acceleration.y = factor1 * total_tidal_force_y;
            particle.tides.parameters.output.acceleration.z = factor1 * total_tidal_force_z;
        }
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
    //particle.tides.parameters.output.acceleration.x += factor2 * sum_total_tidal_force.x;
    //particle.tides.parameters.output.acceleration.y += factor2 * sum_total_tidal_force.y;
    //particle.tides.parameters.output.acceleration.z += factor2 * sum_total_tidal_force.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    // Also factorize
    tidal_host_particle.tides.parameters.output.acceleration.x =
        -1.0 * factor2 * sum_total_tidal_force.x;
    tidal_host_particle.tides.parameters.output.acceleration.y =
        -1.0 * factor2 * sum_total_tidal_force.y;
    tidal_host_particle.tides.parameters.output.acceleration.z =
        -1.0 * factor2 * sum_total_tidal_force.z;
}
