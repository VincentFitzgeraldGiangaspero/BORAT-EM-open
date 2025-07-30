import pickle

def save(filename, *args):
    # Get global dictionary
    glob = globals()
    d = {}
    for v in args:
        # Copy over desired values
        d[v] = glob[v]
    with open(filename, 'wb') as f:
        # Put them in the file 
        pickle.dump(d, f)

def load(filename):
    # Get global dictionary
    glob = globals()
    with open(filename, 'rb') as f:
        for k, v in pickle.load(f).items():
            # Set each global variable to the value from the file
            glob[k] = v
    
def exportSolution(filename,variable):
    with open(filename, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([variable], f)
        f.close

def exportSolutionList(filename,list):
    with open(filename, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump(list, f)
        f.close
        
# Getting back the objects:
def importSolutionList(filename):
    with open(filename,'rb') as f:  # Python 3: open(..., 'rb')
        variable = pickle.load(f)
    return variable

