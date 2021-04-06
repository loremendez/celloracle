# celloracle
testing celloracle using a Pancreas dataset

Originally, the package was written so we propagate the signal of the GRN to infer the direction of movement when causing a perturbation. I modified the library so we can propagate the signal without perturbation. </br>

Se need to modify the *_do_simulation_* function in 'celloracle/trajectory/oracle_GRN.py':

    def _do_simulation(coef_matrix, simulation_input, gem, n_propagation):

        delta_input = simulation_input - gem

        delta_simulated = delta_input.copy()
        for i in range(n_propagation):

            if np.all((delta_simulated == 0)):
                gem_tmp = gem.dot(coef_matrix)
            else:
                delta_simulated = delta_simulated.dot(coef_matrix)
                delta_simulated[delta_input != 0] = delta_input
                gem_tmp = gem + delta_simulated

            gem_tmp[gem_tmp<0] = 0
            delta_simulated = gem_tmp - gem

        gem_simulated = gem + delta_simulated

        return gem_simulated
