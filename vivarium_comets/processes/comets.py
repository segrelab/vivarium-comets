'''
Execute by running: ``python template/processes/template_process.py``

TODO: Replace the template code to implement your own process.
'''

import os
import pytest

import matplotlib.pyplot as plt
import cometspy as c
import cobra.io

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
        'models': [],
        'metabolite_ids': [],
        'defaultVmax': 18.5,
        'defaultKm': 0.000015,
        'maxCycles': 1,
        'spaceWidth': 1,
        'maxSpaceBiomass': 10,
        'minSpaceBiomass': 1e-11,
        'writeMediaLog': True,
        'MediaLogRate': 1,

    }

    def __init__(self, parameters=None):
        super().__init__(parameters)
        dimensions = self.parameters['dimensions']
        models = self.parameters['models']
        self.model_ids = [model.id for model in models]
        metabolite_ids = self.parameters['metabolite_ids']
        defaultVmax = self.parameters['defaultVmax']
        defaultKm = self.parameters['defaultKm']
        maxCycles = self.parameters['maxCycles']
        spaceWidth = self.parameters['spaceWidth']
        maxSpaceBiomass = self.parameters['maxSpaceBiomass']
        minSpaceBiomass = self.parameters['minSpaceBiomass']
        writeMediaLog = self.parameters['writeMediaLog']
        MediaLogRate = self.parameters['MediaLogRate']
        
    def ports_schema(self):
        return {
            # Declare a port for biomass
            'Biomass': {
                # Define the variable that goees through that port
                model_id: {
                    # Define how that variable operates
                    '_default': 0.0,
                    '_updater': 'set', # Use set, since we are not returning the delta but the new value
                    '_emit': True
                } for model_id in self.model_ids
            },
            
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
        biomass = states['Biomass']['e_coli_core']
        metabolites = states['Metabolites']
        
        print(metabolites)
        
        # Run COMETS (need to use timestep somewhere in here)
        ##################################################################
        # Create empty 1x1 layout
        test_tube = c.layout() # TODO: Actually use the dimensions parameter
        
        # Add all the metabolites
        # Kludge- loop through all the metabolites and set specific, only if greater than 0
        for met_id in metabolites:
            if metabolites[met_id] > 0:
                test_tube.set_specific_metabolite(met_id, metabolites[met_id])

        # Set the models' initial biomass and add it to the layout
        for model in self.parameters['models']:
            # First two numbers are the x & y coordinates of the COMETS 2D grid
            # COMETS uses 0 indexing, so 0 0 is the first square
            model.initial_pop = [0, 0, biomass]
            test_tube.add_model(model)
        
        # Hardcode the simulation parameters
        sim_params = c.params()
        sim_params.set_param('defaultVmax', self.parameters['defaultVmax'])
        sim_params.set_param('defaultKm', self.parameters['defaultKm'])
        sim_params.set_param('maxCycles', self.parameters['maxCycles'])
        sim_params.set_param('timeStep', timestep)
        sim_params.set_param('spaceWidth', self.parameters['spaceWidth'])
        sim_params.set_param('maxSpaceBiomass', self.parameters['maxSpaceBiomass'])
        sim_params.set_param('minSpaceBiomass', self.parameters['minSpaceBiomass'])
        sim_params.set_param('writeMediaLog', self.parameters['writeMediaLog'])
        sim_params.set_param('MediaLogRate', self.parameters['MediaLogRate'])

        # Define an experiment
        experiment = c.comets(test_tube, sim_params)
    
        # Run the simulation
        experiment.run()
        ##################################################################
        
        # Get the next biomass and metabolites from the COMETS run
        next_biomass = experiment.total_biomass['e_coli_core'][1] # FIXME: Don't hardocde the model name
        print(next_biomass)

        media = experiment.media.copy()
        print(media)
        previous_media = media.loc[media['cycle'] == media['cycle'].max()-1]
        most_recent_media = media.loc[media['cycle'] == media['cycle'].max()]
        # If a metabolite is not in the most recent media, but was in the previous media, set it to 0
        # This is a kludge to deal with the fact that COMETS doesn't report metabolites that are not present
        if len(most_recent_media) < len(previous_media): # FIXME: This won't always work if you have added and removed metabolites in the same timestep
            # Get the metabolites that are in the previous media but not the most recent media
            missing_metabolites = set(previous_media['metabolite']) - set(most_recent_media['metabolite'])
            # Add rows for the missing metabolites
            for met_id in missing_metabolites:
                most_recent_media = most_recent_media.append({'metabolite': met_id, 'conc_mmol': 0}, ignore_index=True) #TODO: Need to make space explicit
        next_metabolite = most_recent_media.set_index('metabolite').to_dict()['conc_mmol']
        print(next_metabolite)
        
        return {
            'Biomass': {
                'e_coli_core': next_biomass
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
    e_coli_cobra = cobra.io.load_model('textbook')

    # Translate the cobra format into the comets format
    # This is using COMETSpy outside of the process, is that okay?
    e_coli = c.model(e_coli_cobra)

    # remove the bounds from glucose import (will be set dynamically by COMETS)
    # The bounds will be over written by the Michaelis-Menten kinetics
    # By default the bounds are 0 and 1000, can cause problems
    e_coli.change_bounds('EX_glc__D_e', -1000, 1000)
    e_coli.change_bounds('EX_ac_e', -1000, 1000)
    
    # Set the settings
    comets_config = {'time_step': 1.0,
                     'dimensions': [1,1],
                     'models': [e_coli],
                     'metabolite_ids': [met.id for met in e_coli_cobra.metabolites]}
    comets_sim_settings = {
        'experiment_id': 'foo'}

    # Declare the initial state, mirroring the ports structure
    comets_initial_state = {
        'Biomass': {'e_coli_core': 0.000005},
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
    plt.plot(comets_output['time'], comets_output['Biomass']['e_coli_core'])
    plt.savefig('biomass.png')
    
    plt.clf()
    
    plt.plot(comets_output['time'], comets_output['Metabolites']['glc__D_e'], label = "glucose")
    plt.plot(comets_output['time'], comets_output['Metabolites']['ac_e'], label = "acetate")
    plt.plot(comets_output['time'], comets_output['Metabolites']['etoh_e'], label = "ethanol")
    plt.plot(comets_output['time'], comets_output['Metabolites']['for_e'], label = "formate")
    plt.legend()
    plt.savefig('media.png')


# run module with python vivarium_comets/processes/comets.py
if __name__ == '__main__':
    main()
