import os
import inspect
import shutil
import filecmp
import posidonius
from posidonius.tests import common
from test_base import TestBase

class Flattening(TestBase):

    def setUp(self):
        TestBase.setUp(self)
        self.current_filename, ignore = os.path.splitext(os.path.basename(__file__)) # Filename without extension
        self.current_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def tearDown(self):
        TestBase.tearDown(self)

    def test_solar_like(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_tides = True
        consider_disk_interaction = False
        consider_rotational_flattening = True
        #consider_general_relativity = False
        consider_general_relativity = "Kidder1995" # Assumes one central massive body
        #consider_general_relativity = "Anderson1975" # Assumes one central massive body
        #consider_general_relativity = "Newhall1983" # Considers all bodies
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity)

        star_mass = 1.0 # Solar masses
        star_rotation_period = 8. # hours
        star_dissipation_factor_scale = 1.0
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        star_evolution_type = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution_type = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        #star_evolution_type = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        #star_evolution_type = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution_type = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution_type = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution_type = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        #star_evolution_type = posidonius.NonEvolving()
        universe.add_solar_like(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_rotation_period, star_evolution_type, wind_k_factor=4.0e-18, wind_rotation_saturation=1.7592918860102842)
        common.basic_configuration(universe)

        ############################################################################
        whfast_alternative_coordinates="Jacobi"
        #whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        #universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        self.assertTrue(filecmp.cmp(json_filename, expected_json_filename, shallow=False), "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))
        shutil.rmtree(os.path.dirname(json_filename))


    def test_brown_dwarf(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_tides = True
        consider_disk_interaction = False
        consider_rotational_flattening = True
        #consider_general_relativity = False
        consider_general_relativity = "Kidder1995" # Assumes one central massive body
        #consider_general_relativity = "Anderson1975" # Assumes one central massive body
        #consider_general_relativity = "Newhall1983" # Considers all bodies
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity)

        star_mass = 0.08 # Solar masses
        star_dissipation_factor_scale = 1.0
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        #star_evolution_type = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution_type = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        star_evolution_type = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        #star_evolution_type = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution_type = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution_type = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution_type = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        #star_evolution_type = posidonius.NonEvolving()
        universe.add_brown_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_evolution_type)
        common.basic_configuration(universe)

        ############################################################################
        whfast_alternative_coordinates="Jacobi"
        #whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        #universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        self.assertTrue(filecmp.cmp(json_filename, expected_json_filename, shallow=False), "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))
        shutil.rmtree(os.path.dirname(json_filename))


    def test_m_dwarf(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_tides = True
        consider_disk_interaction = False
        consider_rotational_flattening = True
        #consider_general_relativity = False
        consider_general_relativity = "Kidder1995" # Assumes one central massive body
        #consider_general_relativity = "Anderson1975" # Assumes one central massive body
        #consider_general_relativity = "Newhall1983" # Considers all bodies
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_tides, consider_rotational_flattening, consider_disk_interaction, consider_general_relativity)

        star_mass = 0.10 # Solar masses
        star_rotation_period = 70. # hours
        star_dissipation_factor_scale = 1.0
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        star_evolution_type = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution_type = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        #star_evolution_type = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        #star_evolution_type = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution_type = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution_type = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution_type = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        #star_evolution_type = posidonius.NonEvolving()
        universe.add_m_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, star_rotation_period, star_evolution_type)
        common.basic_configuration(universe)

        ############################################################################
        whfast_alternative_coordinates="Jacobi"
        #whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        #universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        self.assertTrue(filecmp.cmp(json_filename, expected_json_filename, shallow=False), "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))
        shutil.rmtree(os.path.dirname(json_filename))


