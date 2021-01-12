from __future__ import print_function
def Adk_aff_chrom_analysis(a, b):
    
    #Input data#
    data = a
    elution_start_fraction = "A4"
    elution_end_fraction = "A10"
    eluent_fractions = ['A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10']
    peak_AU_sensitivity = b #What the user wants the minimum peak absorption highlighted to be: smaller is more sensitive#
    blank_data = "Adk 4 Blank Data Example.csv" #Please do not change unless using different sample loop size. Is the control data
    Protein_Extinction_Coefficient = 10430 #cm^-1 M^-1

    #Reading in sample absorption data and formatting
    import pandas as pd
    data = pd.read_csv(data, skiprows = 2)
    absorption_data = data[["ml", " mAU"]]
    
    #Creating dataframe of fraction volumes
    fractions = data[["(Fractions)", "ml.2"]]
    fractions = fractions.dropna()
    fractions.set_index("(Fractions)", inplace=True)
    collected_fractions = fractions.loc[["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13"]]

    #Pulling the eluent start volume and end volume
    elution_start_vol = collected_fractions.at[elution_start_fraction, "ml.2"]

    collected_fractions.reset_index(inplace = True)
    elution_end_fraction_index = collected_fractions.iloc[collected_fractions[collected_fractions["(Fractions)"] == elution_end_fraction].index + 1]
    elution_end_fraction_index.reset_index(inplace = True)
    elution_end_fraction_index = pd.DataFrame(elution_end_fraction_index)
    elution_end_vol = elution_end_fraction_index.at[0, "ml.2"]

    collected_fractions.set_index("(Fractions)", inplace=True)

    #Calculating the peaks and then limiting these to only the eluted region
    absorption_peaks_full = absorption_data[absorption_data[" mAU"] > peak_AU_sensitivity]
    absorption_peaks_elute = absorption_peaks_full[(absorption_peaks_full["ml"] > elution_start_vol) & (absorption_peaks_full["ml"] < elution_end_vol)]
    
    #Plot a non-normalised plot of absorption with elution peaks
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots()
    absorption_data.plot("ml", " mAU", color="blue", ax=ax)
    absorption_peaks_elute.plot("ml", " mAU", color="red", ax=ax)
    ax.set_ylabel("Absorption at 280 nm (AU)")
    ax.set_xlabel("Elution Volume (ml)")
    ax.get_legend().remove()
    ax.set_title("Absorption Spectrum of Affinity Chromatography")
    plt.show()
    
    #Reading in control data and then filtering for only absorption data
    blank_data = pd.read_csv(blank_data, skiprows = 2)
    blank_absorption_data = blank_data[["ml", " mAU"]]
    
    #Calculating normalised data and filtering to only the eluted fractions
    normalised_data = absorption_data
    normalised_data[" mAU"] = normalised_data[" mAU"] - blank_absorption_data[" mAU"]
    normalised_elution_region = normalised_data[normalised_data["ml"] > elution_start_vol]
    normalised_elution_region = normalised_elution_region[normalised_elution_region["ml"] < elution_end_vol]
    
    #Calculate the normalised peak values and plotting
    normalised_elution_peaks = normalised_elution_region[normalised_elution_region[" mAU"] >= peak_AU_sensitivity]

    fig,ax2 = plt.subplots()
    normalised_elution_region.plot("ml", " mAU", color="blue", ax=ax2)
    normalised_elution_peaks.plot("ml", " mAU", color="red", ax=ax2)
    ax2.set_ylabel("Absorption at 280 nm (AU)")
    ax2.set_xlabel("Elution Volume (ml)")
    ax2.set_title("Normalised Absorption Spectrum of Affinity Chromatography (eluent only)")
    ax2.get_legend().remove()
    plt.show()
    
    #Creating dataframe of only the flowthrough and of fractions AFTER eluent, respectfully
    flowthrough = data[["ml", " mAU"]]
    flowthrough = flowthrough[flowthrough["ml"] <= elution_start_vol]
    after_eluent = normalised_data[["ml", " mAU"]]
    after_eluent = after_eluent[after_eluent["ml"] >= elution_end_vol]
    
    #Concatinating the flowthrough set and normalised eluent set and plotting with peaks
    flowthrough_normalised_eluent = pd.concat([flowthrough, normalised_elution_region, after_eluent])

    fig,ax4 = plt.subplots()
    flowthrough_normalised_eluent.plot("ml", " mAU", color="blue", ax=ax4)
    normalised_elution_peaks.plot("ml", " mAU", color="red", ax=ax4)
    ax4.set_ylabel("Absorption at 280 nm (AU)")
    ax4.set_xlabel("Elution Volume (ml)")
    ax4.set_title("Normalised Absorption Spectrum of Affinity Chromatography")
    ax4.get_legend().remove()
    plt.show()
    
    
    #Calculating the start volume of each fraction
    import numpy as np
    from scipy.integrate import simps
    from numpy import trapz

    A1_start = collected_fractions.at["A1", "ml.2"]
    A2_start = collected_fractions.at["A2", "ml.2"]
    A3_start = collected_fractions.at["A3", "ml.2"]
    A4_start = collected_fractions.at["A4", "ml.2"]
    A5_start = collected_fractions.at["A5", "ml.2"]
    A6_start = collected_fractions.at["A6", "ml.2"]
    A7_start = collected_fractions.at["A7", "ml.2"]
    A8_start = collected_fractions.at["A8", "ml.2"]
    A9_start = collected_fractions.at["A9", "ml.2"]
    A10_start = collected_fractions.at["A10", "ml.2"]
    A11_start = collected_fractions.at["A11", "ml.2"]
    A12_start = collected_fractions.at["A12", "ml.2"]
    A13_start = collected_fractions.at["A13", "ml.2"]

    fractions.reset_index(inplace = True)
    A13_end_index = fractions.iloc[fractions[fractions["(Fractions)"] == "A13"].index + 1]
    A13_end_index.reset_index(inplace = True)
    A13_end_index = pd.DataFrame(A13_end_index)
    A13_end = A13_end_index.at[0, "ml.2"]
    fractions.set_index("(Fractions)", inplace=True)

    #Calculating the fraction volume regions and associated absorption values
    fraction_A1 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A1_start) & (flowthrough_normalised_eluent['ml'] < A2_start)]
    fraction_A2 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A2_start) & (flowthrough_normalised_eluent['ml'] < A3_start)]
    fraction_A3 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A3_start) & (flowthrough_normalised_eluent['ml'] < A4_start)]
    fraction_A4 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A4_start) & (flowthrough_normalised_eluent['ml'] < A5_start)]
    fraction_A5 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A5_start) & (flowthrough_normalised_eluent['ml'] < A6_start)]
    fraction_A6 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A6_start) & (flowthrough_normalised_eluent['ml'] < A7_start)]
    fraction_A7 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A7_start) & (flowthrough_normalised_eluent['ml'] < A8_start)]
    fraction_A8 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A8_start) & (flowthrough_normalised_eluent['ml'] < A9_start)]
    fraction_A9 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A9_start) & (flowthrough_normalised_eluent['ml'] < A10_start)]
    fraction_A10 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A10_start) & (flowthrough_normalised_eluent['ml'] < A11_start)]
    fraction_A11 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A11_start) & (flowthrough_normalised_eluent['ml'] < A12_start)]
    fraction_A12 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A12_start) & (flowthrough_normalised_eluent['ml'] < A13_start)]
    fraction_A13 = flowthrough_normalised_eluent[(flowthrough_normalised_eluent['ml'] >= A13_start) & (flowthrough_normalised_eluent['ml'] < A13_end)]
    
    #Create numpy arrays of absorption values in each fraction
    A1_mAU = np.array(fraction_A1[" mAU"])
    A1_ml = np.array(fraction_A1["ml"])

    A2_mAU = np.array(fraction_A2[" mAU"])
    A2_ml = np.array(fraction_A2["ml"])

    A3_mAU = np.array(fraction_A3[" mAU"])
    A3_ml = np.array(fraction_A3["ml"])

    A4_mAU = np.array(fraction_A4[" mAU"])
    A4_ml = np.array(fraction_A4["ml"])

    A5_mAU = np.array(fraction_A5[" mAU"])
    A5_ml = np.array(fraction_A5["ml"])

    A6_mAU = np.array(fraction_A6[" mAU"])
    A6_ml = np.array(fraction_A6["ml"])

    A7_mAU = np.array(fraction_A7[" mAU"])
    A7_ml = np.array(fraction_A7["ml"])

    A8_mAU = np.array(fraction_A8[" mAU"])
    A8_ml = np.array(fraction_A8["ml"])

    A9_mAU = np.array(fraction_A9[" mAU"])
    A9_ml = np.array(fraction_A9["ml"])

    A10_mAU = np.array(fraction_A10[" mAU"])
    A10_ml = np.array(fraction_A10["ml"])

    A11_mAU = np.array(fraction_A11[" mAU"])
    A11_ml = np.array(fraction_A11["ml"])
    
    A12_mAU = np.array(fraction_A12[" mAU"])
    A12_ml = np.array(fraction_A12["ml"])

    A13_mAU = np.array(fraction_A13[" mAU"])
    A13_ml = np.array(fraction_A13["ml"])
    
    
    #Calculate areas by two different methods
    A1_trapz_area = trapz(A1_mAU, A1_ml)
    A1_simps_area = simps(A1_mAU, A1_ml)
    
    A2_trapz_area = trapz(A2_mAU, A2_ml)
    A2_simps_area = simps(A2_mAU, A2_ml)
    
    A3_trapz_area = trapz(A3_mAU, A3_ml)
    A3_simps_area = simps(A3_mAU, A3_ml)
    
    A4_trapz_area = trapz(A4_mAU, A4_ml)
    A4_simps_area = simps(A4_mAU, A4_ml)
    
    A5_trapz_area = trapz(A5_mAU, A5_ml)
    A5_simps_area = simps(A5_mAU, A5_ml)
    
    A6_trapz_area = trapz(A6_mAU, A6_ml)
    A6_simps_area = simps(A6_mAU, A6_ml)
    
    A7_trapz_area = trapz(A7_mAU, A7_ml)
    A7_simps_area = simps(A7_mAU, A7_ml)
    
    A8_trapz_area = trapz(A8_mAU, A8_ml)
    A8_simps_area = simps(A8_mAU, A8_ml)
    
    A9_trapz_area = trapz(A9_mAU, A9_ml)
    A9_simps_area = simps(A9_mAU, A9_ml)
    
    A10_trapz_area = trapz(A10_mAU, A10_ml)
    A10_simps_area = simps(A10_mAU, A10_ml)
    
    A11_trapz_area = trapz(A11_mAU, A11_ml)
    A11_simps_area = simps(A11_mAU, A11_ml)
    
    A12_trapz_area = trapz(A12_mAU, A12_ml)
    A12_simps_area = simps(A12_mAU, A12_ml)

    A13_trapz_area = trapz(A13_mAU, A13_ml)
    A13_simps_area = simps(A13_mAU, A13_ml)
    
    #Create DataFrame table of results and calculate average areas and fraction sizes
    fraction_areas = [['A1', A1_trapz_area, A1_simps_area, ((A1_trapz_area + A1_simps_area)/2), (A2_start - A1_start)],
                  ['A2', A2_trapz_area, A2_simps_area, ((A2_trapz_area + A2_simps_area)/2), (A3_start - A2_start)],
                  ['A3', A3_trapz_area, A3_simps_area, ((A3_trapz_area + A3_simps_area)/2), (A4_start - A3_start)],
                  ['A4', A4_trapz_area, A4_simps_area, ((A4_trapz_area + A4_simps_area)/2), (A5_start - A4_start)],
                  ['A5', A5_trapz_area, A5_simps_area, ((A5_trapz_area + A5_simps_area)/2), (A6_start - A5_start)],
                  ['A6', A6_trapz_area, A6_simps_area, ((A6_trapz_area + A6_simps_area)/2), (A7_start - A6_start)],
                  ['A7', A7_trapz_area, A7_simps_area, ((A7_trapz_area + A7_simps_area)/2), (A8_start - A7_start)],
                  ['A8', A8_trapz_area, A8_simps_area, ((A8_trapz_area + A8_simps_area)/2), (A9_start - A8_start)],
                  ['A9', A9_trapz_area, A9_simps_area, ((A9_trapz_area + A9_simps_area)/2), (A10_start - A9_start)],
                  ['A10', A10_trapz_area, A10_simps_area, ((A10_trapz_area + A10_simps_area)/2), (A11_start - A10_start)],
                  ['A11', A11_trapz_area, A11_simps_area, ((A11_trapz_area + A11_simps_area)/2), (A12_start - A11_start)],
                  ['A12', A12_trapz_area, A12_simps_area, ((A12_trapz_area + A12_simps_area)/2), (A13_start - A12_start)],
                  ['A13', A13_trapz_area, A13_simps_area, ((A13_trapz_area + A13_simps_area)/2), (A13_end - A13_start)]]

    fraction_areas = pd.DataFrame(fraction_areas, columns = ['Fraction', 'Trapz Area', 'Simps Area', 'Mean Area', 'Fraction Size']) 
    fraction_areas.set_index("Fraction", inplace=True)
    
    #Calculation of protein concentration, appending to output table and plotting bar chart of eluent fraction protein concentrations
    #Beer lambert law is A = elc where e= extinction coefficient, l=pathlength, c=concentration, A=absorption
    #UV flow cell is 2 mm (0.002 cm) in diameter (pathlength). BslA-1 extinction coefficient = 9970 cm^-1 M^-1 calculated at http://www.biomol.net/en/tools/proteinextinction.htm
    
    A1_mean_absorption_per_ml = (fraction_areas.at['A1', 'Mean Area']/1000) / (fraction_areas.at['A1', 'Fraction Size'])
    A1_protein_concentration_mM = (A1_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A2_mean_absorption_per_ml = (fraction_areas.at['A2', 'Mean Area']/1000) / (fraction_areas.at['A2', 'Fraction Size'])
    A2_protein_concentration_mM = (A2_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A3_mean_absorption_per_ml = (fraction_areas.at['A3', 'Mean Area']/1000) / (fraction_areas.at['A3', 'Fraction Size'])
    A3_protein_concentration_mM = (A3_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A4_mean_absorption_per_ml = (fraction_areas.at['A4', 'Mean Area']/1000) / (fraction_areas.at['A4', 'Fraction Size'])
    A4_protein_concentration_mM = (A4_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A5_mean_absorption_per_ml = (fraction_areas.at['A5', 'Mean Area']/1000) / (fraction_areas.at['A5', 'Fraction Size'])
    A5_protein_concentration_mM = (A5_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A6_mean_absorption_per_ml = (fraction_areas.at['A6', 'Mean Area']/1000) / (fraction_areas.at['A6', 'Fraction Size'])
    A6_protein_concentration_mM = (A6_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A7_mean_absorption_per_ml = (fraction_areas.at['A7', 'Mean Area']/1000) / (fraction_areas.at['A7', 'Fraction Size'])
    A7_protein_concentration_mM = (A7_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A8_mean_absorption_per_ml = (fraction_areas.at['A8', 'Mean Area']/1000) / (fraction_areas.at['A8', 'Fraction Size'])
    A8_protein_concentration_mM = (A8_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A9_mean_absorption_per_ml = (fraction_areas.at['A9', 'Mean Area']/1000) / (fraction_areas.at['A9', 'Fraction Size'])
    A9_protein_concentration_mM = (A9_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    

    A10_mean_absorption_per_ml = (fraction_areas.at['A10', 'Mean Area']/1000) / (fraction_areas.at['A10', 'Fraction Size'])
    A10_protein_concentration_mM = (A10_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    
    A11_mean_absorption_per_ml = (fraction_areas.at['A11', 'Mean Area']/1000) / (fraction_areas.at['A11', 'Fraction Size'])
    A11_protein_concentration_mM = (A11_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    
    A12_mean_absorption_per_ml = (fraction_areas.at['A12', 'Mean Area']/1000) / (fraction_areas.at['A12', 'Fraction Size'])
    A12_protein_concentration_mM = (A12_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    
    A13_mean_absorption_per_ml = (fraction_areas.at['A13', 'Mean Area']/1000) / (fraction_areas.at['A13', 'Fraction Size'])
    A13_protein_concentration_mM = (A13_mean_absorption_per_ml / ( Protein_Extinction_Coefficient * 0.002)) *1000
    
    fraction_areas["Mean Protein Concentration (mM)"] = [A1_protein_concentration_mM,
                                                         A2_protein_concentration_mM,
                                                         A3_protein_concentration_mM,
                                                         A4_protein_concentration_mM,
                                                         A5_protein_concentration_mM,
                                                         A6_protein_concentration_mM,
                                                         A7_protein_concentration_mM,
                                                         A8_protein_concentration_mM,
                                                         A9_protein_concentration_mM,
                                                         A10_protein_concentration_mM,
                                                         A11_protein_concentration_mM,
                                                         A12_protein_concentration_mM,
                                                         A13_protein_concentration_mM]
    
    print(fraction_areas)
    eluent_fractions = pd.DataFrame(eluent_fractions, columns = ['Eluent Fractions'])
    eluent_fraction_areas = fraction_areas.loc[eluent_fractions['Eluent Fractions']]
    eluent_fraction_areas.plot(y = 'Mean Protein Concentration (mM)', kind = 'bar', xlabel = "Fraction", ylabel = "Mean Protein Concentration (mM)", title = "Eluent Fraction Protein Concentrations", legend = None)


#Calling the function

Adk_aff_chrom_analysis("Adk1 500mL Culture Data.csv", 350)

Adk_aff_chrom_analysis("Adk4 Arctic Express 500mL Culture Data.csv", 250)