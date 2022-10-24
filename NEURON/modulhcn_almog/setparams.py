from neuron import h
import protocol

defVals = protocol.get_defvals()

def setparams(params):
  params_copy = params.copy()
  keys = params.keys()
  for ikey in range(0,len(keys)):
    params_copy[keys[ikey]] = params[keys[ikey]]*defVals[keys[ikey]]

  print("setparams: params_copy = "+str(params_copy))
  for ikey in range(0,len(keys)):
    key = keys[ikey]
    print(key+" = "+str(params_copy[key]))
    h(key+" = "+str(params_copy[key]))
  h("tfunk()")
