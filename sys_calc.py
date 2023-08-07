import pandas as pd
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

list_of_quantities = [' Gamma112_EPD', ' Gamma132_EPD', ' Gamma112_EPD1', ' Gamma132_EPD1', ' Gamma112_EPD_ESE', 
                        ' Gamma132_EPD_ESE', ' Gamma112_EPD1_ESE', ' Gamma132_EPD1_ESE', ' kappa112_EPD', ' kappa132_EPD', 
                        ' kappa112_EPD1', ' kappa132_EPD1', ' kappa112_EPD_ESE', ' kappa132_EPD_ESE', ' kappa112_EPD1_ESE', 
                        ' kappa132_EPD1_ESE']

def plotGraph(Gamma_Default, final_df, set, quant, R):
    x = ['sys_0', 'sys_1', 'sys_2', 'sys_3', 'sys_4', 'sys_5', 'sys_6', 'sys_7']
    x2 = ['tot_sys', 'sys_1', 'sys_2', 'sys_3', 'sys_4', 'sys_5', 'sys_6', 'sys_7']
    y = []
    y_error = []
    fig = plt.figure(figsize=(12, 8))
    #fig, axs = plt.subplots(nrows=9, , sharex=True)

    #for j in range(0 + 3*set, 3 + 3*set):
        #y_tmp = [Gamma_Default[list_of_quantities[quant] + "_" + str(i)][j] for i in [0,1,2,3,4,5,6] ]
        #y_tmp.append(Gamma_Default[list_of_quantities[quant] + "_0"][j])
        #y.append( y_tmp )

        #y_error_tmp = [Gamma_Default[list_of_quantities[quant] + "_err_" + str(i)][j] for i in [0,1,2,3,4,5,6] ]
        #y_error_tmp.append(final_df[list_of_quantities[quant] + "_sys"][j])
        #y_error.append( y_error_tmp ) #range(0, 8)]
        
        #ax = fig.add_subplot(3, 1, j - 3*set +1)
        #ax.errorbar(x, y[j - 3*set], yerr=y_error[j - 3*set], fmt='o')
        #ax.set_title(list_of_quantities[quant] + ' Centrality ' + str(j))

    cent = 2 * set

    y_tmp = [Gamma_Default[list_of_quantities[quant] + '_' + str(i) + '_recipe' + str(R)][cent] for i in [0,1,2,3,4,5,6,7] ]
    y_error_tmp = [Gamma_Default[list_of_quantities[quant] + "_err_" + str(i) + '_recipe' + str(R)][cent] for i in [0,1,2,3,4,5,6,7] ]
    ax = fig.add_subplot(2, 2, 1)
    ax.errorbar(x, y_tmp, yerr=y_error_tmp, fmt='o')
    ax.set_title(list_of_quantities[quant] + ' Recipe ' + str(R) + ' Centrality ' + str(cent))

    y_tmp2 = [Gamma_Default[list_of_quantities[quant] + "_" + str(i) + '_recipe' + str(R)][cent] for i in [0,1,2,3,4,5,6,7] ]
    y_error_tmp2 = []
    for i in range(0, 8):
        if i == 0:
            y_error_tmp2.append(final_df[list_of_quantities[quant] + '_recipe' + str(R) + '_sys'][cent])
        else:
            y_error_tmp2.append(Gamma_Default[list_of_quantities[quant] + '_recipe' + str(R) + '_sys_' + str(i)][cent])
    
    ax = fig.add_subplot(2, 2, 3)
    ax.errorbar(x2, y_tmp2, yerr=y_error_tmp2, fmt='o')
    ax.set_title(list_of_quantities[quant] + ' Sys. Err. Contribution for Recipe' + str(R) + ' Cent. ' + str(cent))

    cent = 2 * set + 1
    if cent >= 9:
        return fig

    y_tmp = [Gamma_Default[list_of_quantities[quant] + "_" + str(i) + '_recipe' + str(R)][cent] for i in [0,1,2,3,4,5,6,7] ]
    y_error_tmp = [Gamma_Default[list_of_quantities[quant] + "_err_" + str(i) + '_recipe' + str(R)][cent] for i in [0,1,2,3,4,5,6,7] ]
    ax = fig.add_subplot(2, 2, 2)
    ax.errorbar(x, y_tmp, yerr=y_error_tmp, fmt='o')
    ax.set_title(list_of_quantities[quant] + ' Recipe ' + str(R) + ' Centrality ' + str(cent))

    y_tmp2 = [Gamma_Default[list_of_quantities[quant] + "_" + str(i) + '_recipe' + str(R)][cent] for i in [0,1,2,3,4,5,6,7] ]
    y_error_tmp2 = []
    for i in range(0, 8):
        if i == 0:
            y_error_tmp2.append(final_df[list_of_quantities[quant] + '_recipe' + str(R) + '_sys'][cent])
        else:
            y_error_tmp2.append(Gamma_Default[list_of_quantities[quant] + '_recipe' + str(R) + '_sys_' + str(i)][cent])
    
    ax = fig.add_subplot(2, 2, 4)
    ax.errorbar(x2, y_tmp2, yerr=y_error_tmp2, fmt='o')
    ax.set_title(list_of_quantities[quant] + ' Systematic Uncertainty Contribution for Recipe' + str(R) + ' Centrality ' + str(cent))

    return fig


def main_function():

    final_df = pd.DataFrame()
    pp = PdfPages('sys_errs_visualization.pdf')

    for R in range(1, 5):
        Gamma_Default = pd.read_csv('./KFParticle_Results/sys_0/gamma_results_recipe' + str(R) + '.txt')
        Gamma_Default = Gamma_Default.add_suffix('_0_recipe' + str(R))

        for i in range (1, 8):
            Gamma_sys = pd.read_csv('./KFParticle_Results/sys_' + str(i) + '/gamma_results_recipe' + str(R) + '.txt')
            Gamma_Default = Gamma_Default.merge(Gamma_sys.add_suffix('_' + str(i) + '_recipe' + str(R)), left_on='Centrality_0_recipe' + str(R), right_on='Centrality_' + str(i) + '_recipe' + str(R))

            for quant in list_of_quantities:
                #Gamma_Default[quant + '_sys_' + str(i)] = Gamma_Default.apply(lambda x: 
                                                                                    #math.sqrt( (x[quant + '_0'] - x[quant + '_' + str(i)]) ** 2 - 
                                                                                                #abs(x[quant + '_err_0'] ** 2 - x[quant + '_err_' + str(i)] ** 2) ) / math.sqrt(12)  
                                                                                    #if ((x[quant + '_0'] - x[quant + '_' + str(i)]) ** 2 - 
                                                                                                #abs(x[quant + '_err_0'] ** 2 - x[quant + '_err_' + str(i)] ** 2)) > 0 
                                                                                    #else 0, axis=1)
                Gamma_Default[quant + '_recipe' + str(R) + '_sys_' + str(i)] = Gamma_Default.apply(lambda x: (x[quant + '_0_recipe' + str(R)] - x[quant + '_' + str(i) + '_recipe' + str(R)]) / math.sqrt(12), axis=1)

        
        for quant in list_of_quantities:
            final_df[quant + '_recipe' + str(R) + '_sys'] = Gamma_Default.apply(lambda x: math.sqrt( x[quant + '_recipe' + str(R) + '_sys_1'] ** 2 + x[quant + '_recipe' + str(R) + '_sys_2'] ** 2
                                                                                + x[quant + '_recipe' + str(R) + '_sys_3'] ** 2 + x[quant + '_recipe' + str(R) + '_sys_4'] ** 2 
                                                                                + x[quant + '_recipe' + str(R) + '_sys_5'] ** 2 + x[quant + '_recipe' + str(R) + '_sys_6'] ** 2
                                                                                + x[quant + '_recipe' + str(R) + '_sys_7'] ** 2), axis = 1)

    
        # Plotting
    
        for k in range(0, len(list_of_quantities)):
            for l in range(0,5):
                Figure1 = plotGraph(Gamma_Default, final_df, l, k, R)
                pp.savefig(Figure1, orientation='landscape')
                plt.close(Figure1)
    
    pp.close()
    
    print("done")
    
    final_df.to_csv('./KFParticle_Results/sys_errs.csv')

if __name__ == "__main__":
    main_function()