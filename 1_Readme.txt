any Qs email me: collinson.lp@virginmedia.com

Files are alphabetical (sad) but should be clustered by function
the essentials
- hex2rgb and rgb2hex (for PIPs)
- The data for passive pathogen and toxic pathogen (figure 15 and 18). This contains the results from all competitions (every gamma value for every spectrum),  GamESSs for each value of spectrum, and the kipESS. NOTE - for the majority of my data, gamma is scanned between 0.0 and 0.2 with 99 steps, as the ESS almost always lies in this range and this increases the resolution of the ESS calc.
- patch_mods, the system of equations used. NutrSteadyState is the settings, so also needs to be saved. Both of these for the regular model, probio (1 nutrient) and probio (2 nutrients) have different patch_mods and steadystates 

anything with patch_mod is the system of equations which will be used, so needs to be saved. Nutr_steadystate is the settings, so also needs to be saved
the regular model, probio (1 nutrient) and probio (2 nutrients) have different patch_mods and steadystates 

the Code to run is split into separate files, numbered depending on their order
most code is 'ot_3stp' (one toxin, three strains). frequency model is prefaced by 'ot_3stp_freq' and probiotic tests by 'probiotic'

1. setting constants - what it says on tin

2. Prod_scan - running the competitions with every value of gamma for each spectrum value. 'prod_scan_simp' only changes spectrum (used for most figures other than parameter sweeps)
A blank 'prod_scan_parameter' is included with the general template for parameter sweeps.

3. Prod_ESS - calculates the ESS production rate for every value of spectrum. Finds the diagonal (res|res) value with the greatest number of unsuccessful invasions. For the frequency model, Prod_ESS_new calculates the ESS as the point between the fewest invasions 'above' and 'below', as there can be regions of mutual invasibility.

4. Prod_PIP is used to generate PIPs for production rate evolution

5. Spec_scan - runs comps for all values of spectrum (using the calcuated ESSs for the resident spectrum). In the freq model, this can take a while as gamma ESS changes over frequency so all frequency values are run also

6. spec_ESS - calculates the ESS value of kip (affinity for pathogen) by finding the diagnoal (res|res) value with the most unsuccessful invasions

7. Spec_PIP - wow a spectrum PIP

The code for the probiotic does not have this many stages, so is contained in:
setting constants
probio with one nutrient
probio with two nutrients
