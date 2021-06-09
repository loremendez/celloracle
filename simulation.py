import pandas as pd
import numpy as np
import scvelo as scv
from scvelo.tools.transition_matrix import transition_matrix

def signal_prop(oracle, vkey='velocity', present='imputed_count', cluster_col='clusters', n_propagation=1):
    cluster_info = oracle.adata.obs[cluster_col]

    TM = transition_matrix(oracle.adata, vkey=vkey)
    TM = pd.DataFrame(TM.toarray(), columns = oracle.adata.obs.index.values, index=oracle.adata.obs.index.values)

    gem_imputed = oracle.adata.to_df(present)
    gem_future = TM.dot(gem_imputed)
    #adata.layers['gem_future'] = TM.dot(gem_imputed)

    #target genes
    target_genes={}
    for cluster in np.unique(cluster_info):
        target_genes[cluster] = set(oracle.cluster_specific_TFdict[cluster].keys()).intersection(set(oracle.adata.var_names))

    simulation_input = gem_future.copy()
    #just use TFs for the simulation
    simulation_input.loc[:,~gem_future.columns.isin(oracle.active_regulatory_genes)] = gem_imputed.loc[:,~gem_future.columns.isin(oracle.active_regulatory_genes)]

    #simulation
    simulated=[]
    for cluster in np.unique(cluster_info):
        simulation_input_ = simulation_input[cluster_info == cluster] #future value of transcription factors
        gem_ = gem_imputed[cluster_info == cluster] #expression matrix
        #####
        simulation_input.loc[:,~simulation_input.columns.isin(target_genes[cluster])] = gem_imputed.loc[:,~simulation_input.columns.isin(target_genes[cluster])]
        #####
        coef_matrix = oracle.coef_matrix_per_cluster[cluster] #coefficient matrix
        delta_input = simulation_input_ - gem_ #difference in TF expression (from t1 to t0), the other genes are zero
        delta_simulated = delta_input.copy()
        for i in range(n_propagation):
            delta_simulated = delta_simulated.dot(coef_matrix) #propagate the difference in expression of TFs using the coefficients
            #delta_simulated[delta_input != 0] = delta_input #do not predict transcription factor expression
            # gene expression cannot be negative. adjust delta values to make sure that gene expression are not negative values.
            gem_tmp = gem_ + delta_simulated #add the change in expression of other genes, to obtain their future value
            gem_tmp[gem_tmp<0] = 0 #do not allow negative counts
            delta_simulated = gem_tmp - gem_ #new difference in expression
        gem_simulated = gem_ + delta_simulated #add the difference in expression to obtain the future value
        simulated_in_the_cluster = gem_simulated
        simulated.append(simulated_in_the_cluster)
    gem_simulated = pd.concat(simulated, axis=0)
    gem_simulated = gem_simulated.reindex(gem_imputed.index)

    return gem_simulated - gem_imputed
