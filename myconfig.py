'''
List of variables to be adjusted
'''

sequence = 'dGdCdATTmCTAATAGmCdAdGdC'  # Enter  a sequence
#dCdUdATTTGGATGTmCdAdGdC  - Dani
#dGdCdATTmCTAATAGmCdAdGdC -Malat
neutral_losses = ['H2O', '2H', 'NH3', 'CH3']
# Define variables
end_3 = 'OH'
end_5 = 'OH'
ring = 'C5H5O' #Constant
phos = 'HPO3S'  # This means all groups are POH, for internal fragments u added this extra, need to take care when chaning
RNA = 'OH'
DNA = 'H'
modification = 'OC2H2' #change this, OC2H2,OC2H5, OC
min_mass  = 98
max_mass = 3000
max_charges = 9

# Select the sample
sample = DNA

# Enter y/n for the following options

# neutral losses
calc_neutral ='n'
#base loss for precursors are combined with neutral loss therefore need to have them both defined to view them
calc_fragments = 'y'
calc_base = 'n'
calc_internal = 'n'



