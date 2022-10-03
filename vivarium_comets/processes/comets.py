'''
Execute by running: ``python template/processes/template_process.py``

TODO: Replace the template code to implement your own process.
'''

import os

from vivarium.core.process import Process
from vivarium.core.composition import (
    simulate_process,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output


# a global NAME is used for the output directory and for the process name
NAME = 'template'


class Comets(Process):
    
    defaults = {
        'dimensions': [1,1]
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
            
            # Use dictionary comprehension to declare schema for all the metabolites
            'Metabolites': {
#                 'metabolite_dictionary': {
                '*': {
                    '_default': 0.0,
                    '_updater': 'set',
                    '_emit': True
                }
            }
        }
    
    def next_update(self, timestep, states):
        # Parse variables from states
        biomass = states['Biomass']['ecoli_biomass']
        metabolites = states['Metabolites']
        
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

        # Define an experiment
        experiment = c.comets(test_tube, sim_params)
    
        # Run the simulation
        experiment.run()
        ##################################################################
        
        # Get the next biomass and metabolites from the COMETS run
        next_biomass = experiment.total_biomass['e_coli_core'][1]
        next_metabolite = experiment.media.copy()
        print(next_metabolite)
        
        return {
            'Biomass': {
                'ecoli_biomass': next_biomass
            },
            'Metabolites': next_metabolite
        }


# functions to configure and run the process
def run_template_process():
    '''Run a simulation of the process.

    Returns:
        The simulation output.
    '''

    # initialize the process by passing in parameters
    parameters = {}
    template_process = Template(parameters)

    # declare the initial state, mirroring the ports structure
    initial_state = {
        'internal': {
            'A': 0.0
        },
        'external': {
            'A': 1.0
        },
    }

    # run the simulation
    sim_settings = {
        'total_time': 10,
        'initial_state': initial_state}
    output = simulate_process(template_process, sim_settings)

    return output


def test_template_process():
    '''Test that the process runs correctly.

    This will be executed by pytest.
    '''
    output = run_template_process()
    # TODO: Add assert statements to ensure correct performance.


def main():
    '''Simulate the process and plot results.'''
    # make an output directory to save plots
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    output = run_template_process()

    # plot the simulation output
    plot_settings = {}
    plot_simulation_output(output, plot_settings, out_dir)


# run module with python template/processes/template_process.py
if __name__ == '__main__':
    main()
