'''
Execute by running: ``python template/processes/template_process.py``

TODO: Replace the template code to implement your own process.
'''

import os
import pytest

import matplotlib.pyplot as plt
import cometspy as c
import cobra.test

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR,
)

# Process, Deriver, and Composer base classes
from vivarium.core.process import Process, Deriver
from vivarium.core.composer import Composer
from vivarium.core.registry import process_registry

from vivarium.core.engine import Engine, pf
from vivarium.library.units import units


# a global NAME is used for the output directory and for the process name
NAME = 'comets'


class Comets(Process):
    
    defaults = {
        'dimensions': [1,1],
        'metabolite_ids': []
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        dimensions = self.parameters['dimensions']
        
    def ports_schema(self):
        return {
            # Declare a port for biomass
            'Biomass': {
                # Define the variable that goees through that port
                'ecoli_biomass': {
                    # Define how that variable operates
                    '_default': 0.0,
                    '_updater': 'set', # Use set, since we are not returning the delta but the new value
                    '_emit': True}},
            
            # Use dictionary comprehension to declare schema for all the metabolites listed in the initial state
            'Metabolites': {
                mol_id: {
                    '_default': 0.0,
                    '_updater': 'set',
                    '_emit': True
                } for mol_id in self.parameters['metabolite_ids']
            }
        }
    
    def next_update(self, timestep, states):
        # Parse variables from states
        biomass = states['Biomass']['ecoli_biomass']
        metabolites = states['Metabolites']
        
        print(metabolites)
        
        # Run COMETS (need to use timestep somewhere in here)
        ##################################################################
        # Create empty 1x1 layout
        test_tube = c.layout()
        
        # Add all the metabolites
        # Kludge- loop through all the metabolites and set specific
        for met_id in metabolites:
            test_tube.set_specific_metabolite(met_id, metabolites[met_id])
        
        # Hard code loading the E coli model
        e_coli_cobra = cobra.test.create_test_model('textbook')
        # Translate the cobra format into the comets format
        e_coli = c.model(e_coli_cobra)

        # remove the bounds from glucose import (will be set dynamically by COMETS)
        # The bounds will be over written by the Michaelis-Menten kinetics
        # By default the bounds are 0 and 1000, can cause problems
        e_coli.change_bounds('EX_glc__D_e', -1000, 1000)
        e_coli.change_bounds('EX_ac_e', -1000, 1000)
        
        # set the model's initial biomass
        # First two numbers are the x & y coordinates of the COMETS 2D grid
        # COMETS uses 0 indexing, so 0 0 is the first square
        e_coli.initial_pop = [0, 0, biomass]

        # add it to the test_tube
        test_tube.add_model(e_coli)
        
        # Hardcode the simulation parameters
        sim_params = c.params()
        sim_params.set_param('defaultVmax', 18.5)
        sim_params.set_param('defaultKm', 0.000015)
        sim_params.set_param('maxCycles', 1)
        sim_params.set_param('timeStep', timestep)
        sim_params.set_param('spaceWidth', 1)
        sim_params.set_param('maxSpaceBiomass', 10)
        sim_params.set_param('minSpaceBiomass', 1e-11)
        sim_params.set_param('writeMediaLog', True)
        sim_params.set_param('MediaLogRate', 1)

        # Define an experiment
        experiment = c.comets(test_tube, sim_params)
    
        # Run the simulation
        experiment.run()
        ##################################################################
        
        # Get the next biomass and metabolites from the COMETS run
        next_biomass = experiment.total_biomass['e_coli_core'][1]

        media = experiment.media.copy()
        print(media)
        most_recent_media = media.loc[media['cycle'] == media['cycle'].max()]
        next_metabolite = most_recent_media.set_index('metabolite').to_dict()['conc_mmol']
        print(next_metabolite)
        
        return {
            'Biomass': {
                'ecoli_biomass': next_biomass
            },
            'Metabolites': next_metabolite
        }


# functions to configure and run the process
def run_comets_process():
    '''Run a the COMETS wrapper for 10 cycles (with timestep of 1.0) with just the E. coli model.

    Returns:
        The simulation output.
    '''
    # Read in the cobra model and get the list of all the metabolites
    e_coli_cobra = cobra.test.create_test_model('textbook')
    
    # Set the settings
    comets_config = {'time_step': 1.0,
                     'dimensions': [1,1],
                     'metabolite_ids': [met.id for met in e_coli_cobra.metabolites]}
    comets_sim_settings = {
        'experiment_id': 'foo'}

    # Declare the initial state, mirroring the ports structure
    comets_initial_state = {
        'Biomass': {'ecoli_biomass': 0.000005},
        'Metabolites': {'glc__D_e': 0.011,
                        'o2_e': 1000,
                        'nh4_e': 1000,
                        'pi_e': 1000,
                        'h2o_e': 1000,
                        'h_e': 1000
                    }
    }

    # Initialize the process
    comets_process = Comets(comets_config)

    # Make the experiment
    comets_exp = Engine(processes={'comets': comets_process},
                    topology={
                        'comets': {'Biomass': ('Biomass',),
                                   'Metabolites': ('Metabolites',)}},
                   initial_state = comets_initial_state)

    # Run the simulation
    comets_exp.update(10.0)

    return comets_exp


def test_comets_process():
    '''Test that the COMETS wrapper process runs correctly.

    This will be executed by pytest.
    '''
    comets_exp = run_comets_process()

    # Check that the experiment has the expected results
    pass


def main():
    # TODO: Update
    '''Simulate the process and plot results.'''
    # Make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    comets_exp = run_comets_process()

    # retrieve the data as a timeseries
    comets_output = comets_exp.emitter.get_timeseries()

    # Plot the simulation output
    plt.plot(comets_output['time'], comets_output['Biomass']['ecoli_biomass'])
    plt.savefig('biomass.png')
    
    plt.clf()
    
    plt.plot(comets_output['time'], comets_output['Metabolites']['glc__D_e'], label = "glucose")
    plt.legend()
    plt.savefig('media.png')


# run module with python vivarium_comets/processes/comets.py
if __name__ == '__main__':
    main()
