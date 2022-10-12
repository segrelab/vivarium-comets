'''
Execute by running: ``python template/processes/template_process.py``

TODO: Replace the template code to implement your own process.
'''

from functools import total_ordering
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
        # TODO: Remove any of the parameters that are not needed to start the simulation
        # Do I need any more parameters? Should I do all of the parameters available in COMETSpy?
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

        # TODO: Start the COMETS simulation so that next update only needs to run one cycle
        
    def ports_schema(self):
        return {
            # Declare a port for biomass
            'Biomass': {
                # Define the variable that goees through that port
                'biomass_grid': {
                    # Define how that variable operates
                    '_default': [],
                    '_updater': 'set', # Use set, since we are not returning the delta but the new value
                    '_emit': True
                }
            },
            
            # Use dictionary comprehension to declare schema for all the metabolites listed in the initial state
            'Metabolites': {
                'metabolite_grid': {
                    '_default': [],
                    '_updater': 'set',
                    '_emit': True
                }
            }
        }
    
    def next_update(self, timestep, states):
        # Parse variables from states
        biomass = states['Biomass']['biomass_grid']
        metabolites = states['Metabolites']['metabolite_grid']
        
        # Run COMETS (need to use timestep somewhere in here)
        ##################################################################
        # Create empty layout from dimensions
        layout = c.layout()
        layout.grid = self.parameters.dimensions

        # Loop through the metabolite grid and set the metabolite for each cell
        for i in range(len(metabolites)):
            # Flip the row index so that the bottom left is (0,0)
            row = len(metabolites) - i - 1
            for j in range(len(metabolites[i])):
                # Loop through the metabolites in the cell
                for met_id in metabolites[i][j]:
                    # If the metabolite is greater than 0, set it
                    if metabolites[i][j][met_id] > 0:
                        layout.set_specific_metabolite_at_location(met_id, [i, row], metabolites[i][j][met_id])


        # Set the models' initial biomass and add it to the layout
        for model in self.parameters['models']:
            # Loop thorugh the biomass grid and set the initial population for each grid cell
            for i, row in enumerate(biomass):
                # Flip the row index so that the bottom left is (0,0)
                row_index = len(biomass) - i - 1
                for j, cell in enumerate(row):
                    # Set the initial biomass only if it is greater than 0
                    if cell[model.id] > 0:
                        model.initial_pop = (row_index, j, cell[model.id])
            # Add the model to the layout
            layout.add_model(model)
        
        # Hardcode the simulation parameters
        sim_params = c.params()
        # TODO: Make this a loop?
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
        experiment = c.comets(layout, sim_params)
    
        # Run the simulation
        experiment.run()
        ##################################################################
        
        # Get the next biomass and metabolites from the COMETS run
        total_biomass = experiment.total_biomass
        current_biomass = total_biomass.loc[total_biomass['cycle'] == total_biomass['cycle'].max()].drop(['cycle'], axis=1).reset_index(drop=True)
        next_biomass = {}
        for col in current_biomass:
            next_biomass[col] = current_biomass[col][0] # This assumes that there is only one row in the dataframe
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
            'Biomass': next_biomass,
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


# functions to configure and run the process
def run_comets_process_2D():
    '''Run a the COMETS wrapper for 10 cycles (with timestep of 1.0) with just the E. coli model in a heterogenous grid.

    Returns:
        The simulation output.
    '''
    # Load in 3 copies of the E. coli model
    model_ecoli1 = cobra.io.read_sbml_model('e_coli_core1.xml')
    model_ecoli2 = cobra.io.read_sbml_model('e_coli_core2.xml')

    # Create COMETS models, set the uptake rate for one of the models to be lower
    ecoli1 = c.model(model_ecoli1)
    ecoli1.id='ecoli1'
    ecoli2 = c.model(model_ecoli2)
    ecoli2.id='ecoli2'
    ecoli2.change_vmax('EX_ac_e',0.0001)

    # Make sure the lower bounds of the uptakes are not zero
    ecoli1.change_bounds('EX_glc__D_e', -1000, 1000)
    ecoli1.change_bounds('EX_ac_e', -1000, 1000)

    ecoli2.change_bounds('EX_glc__D_e', -1000, 1000)
    ecoli2.change_bounds('EX_ac_e', -1000, 1000)
    
    # Set the settings
    comets_config = {'time_step': 0.03,
                     'dimensions': [30,30],
                     'models': [ecoli1, ecoli2],
                     'metabolite_ids': [met.id for met in model_ecoli1.metabolites]} # This only works because all the models have the same metabolites

    # Set up an empty grid for the biomass with a dictionary of 0s in each
    biomass_in_one_cell = {'ecoli1': 0.0, 'ecoli2': 0.0}
    initial_biomass = [[biomass_in_one_cell for i in range(comets_config['dimensions'][0])] for j in range(comets_config['dimensions'][1])]

    # Fill in the biomass grid with the initial biomass
    initial_biomass[15][25]['ecoli1'] = 1e-6
    initial_biomass[15][4]['ecoli2'] = 1e-6

    # Make a dictionary of the initial metabolites
    met_in_one_cell = {'glc__D_e': 0.00005,
                        'o2_e': 1000,
                        'nh4_e': 1000,
                        'pi_e': 1000,
                        'h2o_e': 1000,
                        'h_e': 1000
                    }
    # Make a grid of the initial metabolites
    initial_metabolites = [[met_in_one_cell for i in range(comets_config['dimensions'][0])] for j in range(comets_config['dimensions'][1])]

    # Declare the initial state, mirroring the ports structure
    comets_initial_state = {
        'Biomass': {'biomass_grid': initial_biomass},
        'Metabolites': {'metabolite_grid': initial_metabolites}
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
    comets_exp.update(30.0)

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
