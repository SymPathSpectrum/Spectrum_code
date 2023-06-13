# Spectrum_code
any Qs, email me: collinson.lp@virginmedia.com
using MATLAB version 9.13.0.2049777 (R2022b)

The essentials 
- hex2rgb and rgb2hex (needed for PIPs)
- 3.2.2 and 3.2.3 contains the final abundances, gamma ESSs and Kip ESS for figures 15 and 18
- the patch_mod and NutrSteadyState codes encode the system of equations and steady state function respecitvely, so need to be saved before running comps

--
# Figures 13, 15, 18, S4 and S5 generated using code starting with 'a'
abc_patch_mod and abc_NutrSteadyState need to be saved
1. setting_constants sets the parameters - main difference between the figures
    Fig 13, default, GamP = 0
    Fig 15, same as 13
    Fig 18, default, GamP = 0.2, EP = 2
    
2. a1_production_rate_ESS_and_PIP 
    Fig 13 can be generated here
    13, 15 and 18 all use GamI0 = 0.0, GamIF = 0.2, stpG = 99, kip0 = 3, kipF = 0.05, stpK = 49
    
3. a2_spectrum_ESS_and_PIP
    fig 15 and 18 generated here, no changes needed 
    
-- 
# Figures 11 and S1 also generated using 'a' code
abc_patch_mod and abc_NutrSteadyState need to be saved
1. Setting_constants

2. a3_production_outcomes
    set spectrum as required
    
-- 
# figures S2 and S3 (parameter sweeps) generated using code starting with 'b'
abc_patch_mod and abc_NutrSteadyState need to be saved
1. setting_constants
    S2 and S3, default, GamP = 0
    
2. b1_production_rate_ESS_and_PIP 
    S2 can be generated here, GamI0 = 0.0, GamIF = 0.99, stpG = 99, kip0 = 3, kipF = 0.05, stpK = 1
    S3 needs GamI0 = 0.0, GamIF = 0.99, stpG = 99, kip0 = 3, kipF = 0.05, stpK = 49
    
3. b2_spectrum_ESS_and_PIP
    S3 generated here (one parameter at a time)
    
--
# figures 20, 21 and 22 (frequency) generated with code starting with 'c'
abc_patch_mod and abc_NutrSteadyState need to be saved
1. setting_constants
    Fig 20, default, gamP = 0
    Fig 21 and 22, default, gamP = 0.2, EP = 2

2. c1_production_rate_ESS_and_PIP 
    Figs 20, 21 and 22 generated here
    20, 21 and 22 all use GamI0 = 0.0, GamIF = 0.2, stpG = 99, kip0 = 3, kipF = 0.05, stpK = 49
    
3. c2_spectrum_ESS_and_PIP 
    unrefined as i didn't really get this far!
    Know what you're production rates are doing before trying this step
    
--
# Figure 27 (probiotic, one nutrient) generated using code starting with 'd'
d_patch_mod and d_NutrSteadyState need to be saved
1. d_probiotic_one_nutrient
    fig 27 values generated in heatmaps which reflect the figure
    when looking at 'ideal_value' matrices to get optimal production rates and affinties, be aware that the data is structured differently!
    
--
# Figure 29 (probiotic, two nutrient) generated using code starting with 'e'
e_patch_mod and e_NutrSteadyState need to be saved
1. e_probiotic_two_nutrients
    fig 29 values generated in heatmaps which reflect the figure 
    
--
# Dynamics plots generated using 'indiv_fight' for relevant equation systems ('abc', 'd', or 'e')
