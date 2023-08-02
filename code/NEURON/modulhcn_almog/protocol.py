from neuron import h


defVals = {
  'transvec.x(4)': 118.437,      #gih_end
  'transvec.x(5)': 351.814,      #gih_x2
  'transvec.x(6)': -0.0137979,   #gih_alpha
  'transvec.x(7)': 2.5117,       #gih_start
  'transvec.x(9)': 3.78767,      #gkslow_start
  'transvec.x(10)': -0.0915783,  #gkslow_alpha
  'transvec.x(11)': 205.886,     #gkslow_beta
  'transvec.x(12)': 28.2824,     #gka_start
  'transvec.x(13)': -0.0117721,  #gka_alpha
  'transvec.x(14)': 331.65,      #gka_beta
  'transvec.x(15)': 352.142,     #gna_soma
  'transvec.x(16)': 56.4519,     #gna_api
  'transvec.x(20)': 9.27521e-05, #pcah_soma
  'transvec.x(21)': 0.000155847, #pcah_api
  'transvec.x(25)': 0.00314901,  #pcar_soma
  'transvec.x(26)': 0.000488401, #pcar_api
  'transvec.x(31)': 3.18076,     #gsk_soma
  'transvec.x(32)': 0.524016,    #gsk_dend
  'transvec.x(34)': 0.638741,    #gbk_soma
  'transvec.x(35)': 1.22971,     #gbk_dend
}

variables = [['transvec.x(4)', 0.8, 1.25],
             ['transvec.x(5)', 0.8, 1.25],
             ['transvec.x(6)', 0.8, 1.25],
             ['transvec.x(7)', 0.8, 1.25],
             ['transvec.x(9)', 0.8, 1.25],
             ['transvec.x(10)', 0.8, 1.25],
             ['transvec.x(11)', 0.8, 1.25],
             ['transvec.x(12)', 0.8, 1.25],
             ['transvec.x(13)', 0.8, 1.25],
             ['transvec.x(14)', 0.8, 1.25],
             ['transvec.x(15)', 0.8, 1.25],
             ['transvec.x(16)', 0.8, 1.25],
             ['transvec.x(20)', 0.8, 1.25],
             ['transvec.x(21)', 0.8, 1.25],
             ['transvec.x(25)', 0.8, 1.25],
             ['transvec.x(26)', 0.8, 1.25],
             ['transvec.x(31)', 0.8, 1.25],
             ['transvec.x(32)', 0.8, 1.25],
             ['transvec.x(34)', 0.8, 1.25],
             ['transvec.x(35)', 0.8, 1.25]]



#                   NAME    TYPE      WHERE    AT A FIXED POS?   AMPLITUDE NAME
#                                              OR AT A DISTANCE                
#                                              FROM SOMA?                      
stimulus_types = [ ["st1",  "IClamp", "a_soma", ["fixed", 0.5],    "amp" ],
                   ["syn1", "epsp",   "apic", ["distance", 620], "imax" ] ]

#                   WHAT TYPE OF OUTPUT? VOLTAGE/[CA] AT A
#                   FIXED TIME (LAST TVEC BEFORE GIVEN TIME)
#                   OR MAXIMUM/NSPIKES DURING A GIVEN INTERVAL?
data_storage_types = [ ["fixed", 13000],
                       ["max", [10000,10200] ],
                       ["trace", [9950+5*x for x in range(0,51)] ],
                       ["highrestrace", [9950+x for x in range(0,251)] ],
                       ["highrestraceandspikes", [9950+x for x in range(0,251)] ],
                       ["nspikes", [12000, 15000] ],
                       ["nspikesandothers", [12000, 15000] ] ]


#            STIMULUS_TYPE   AMPLITUDE NAME  PARAMETERS        
#                                                              
#                                                              
stimuli = [ [ 0, [ ["del", 10000], ["dur", 3000] ] ],
            [ 0, [ ["del", 10000], ["dur", 100] ] ],
            [ 1, [ ["onset", 10000], ["tau0", 0.5], ["tau1", 5.0] ] ],
            [ 0, [ ["del", 10000], ["dur", 5000] ] ],
            [ 0, [ ["del", 10000], ["dur", 5] ] ] ]

#                NAME     WHAT    WHERE
#                                   WHICH    AT WHICH
#                                   BRANCH   LOCATIONS
recordings = [ [ "vsoma",  "v",   [ ["a_soma", [0.5] ] ] ],
               [ "vdend",  "v",   [ ["apic", [0.05*x for x in range(0,21)] ], 
                                    ["dend", [0.05*x for x in range(0,21)] ] ] ],
               [ "cadend", "cai", [ ["apic", [0.05*x for x in range(0,21)] ] ] ] ]

#                STIMULUS INDEX
#                        AMPLITUDES
#                                                 RECORDING INDEX
#                                                         DATA STORAGE INDEX
setup =      [ [ [1],   [0.25, 0.5],              0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO A 100 ms DC,
               [ [4],   [0.7],                    0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO A SHORT SOMATIC STIMULUS
               [ [2],   [1.75],                   0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO A SHORT SOMATIC STIMULUS
               [ [4,2], [[0.7, 1.75]],            0,      4 ],            # MEMBRANE POTENTIAL TRACE RESPONSE TO COMBINATION OF SHORT SOMATIC AND APICAL STIMULI
               [ [3],   [0.8, 0.9, 1.0],          0,      6 ] ]           # NUMBERS OF SPIKES, CVISI AND NSPIKESPERBURST AS A RESPONSE TO A LONG DC

stop_if_no_spikes = [ [0, 1], [0], [0], [1], [1, 1, 1] ]
stop_if_more_spikes_or_as_many_as = [ [2, 5], [5], [5], [5], [100, 200, 300] ]  # when nSpikes >= x)

# Return the full matrix of the step-wise fitting variables
def get_variables():
  global variables
  return variables
def get_defvals():
  global defVals
  return defVals

# Return the stimulus array
def get_stimuli():
  global stimuli
  return stimuli

# Return the stimulus types array
def get_stimulus_types():
  global stimulus_types
  return stimulus_types

# Return the data storage array
def get_data_storage_types():
  global data_storage_types
  return data_storage_types

# Return the recordings array
def get_recordings():
  global recordings
  return recordings

# Return the stimulus-recording setup array
def get_setup():
  global setup
  return setup

def get_nspike_restrictions():
  global stop_if_no_spikes, stop_if_more_spikes_or_as_many_as
  return [stop_if_no_spikes, stop_if_more_spikes_or_as_many_as]




