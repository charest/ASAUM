"""
This is a sample plotting utility.

Right now it only plots csv files, but it could be used to plot other files.
"""


################################################################################
# Imports
################################################################################
import argparse
import cycler as cyc
import matplotlib.pyplot as plt
import pandas as pd
import sys


################################################################################
# Declares
################################################################################
def comments(filename):
    cmts = []
    with open(filename) as f:
        for i, line in enumerate(f.readlines()):
            l = line.lstrip()
            if l:
                if l[0] == '#':
                    cmts.append(i)
    return cmts

def labels(s):
    try:
        splt = s.split(',')
        if len(splt)==1:
            return splt[0], splt[0]
        var, label = map(str, splt)
        return var, label
    except:
        raise argparse.ArgumentTypeError("Variables must be of the form: 'varname,label' or 'varname'")


################################################################################
# Main
################################################################################

# setup command line arguments
parser = argparse.ArgumentParser(description='Plot a csv file.')
parser.add_argument('--files', nargs='+', help='List of files to print.')
parser.add_argument('--plotvar', nargs='*', dest='variables', action='append', type=labels, help="Variables to plot.")
parser.add_argument('--xvar', dest='xvar', type=labels, help="Independant variable.")

args = parser.parse_args()

# load all the data files
frames = []
variables = {}
xvar = str()

for f in args.files:
    print('Reading {}'.format(f))
    cmts = comments(f)
    df = pd.read_csv(f, sep='\s*,\s*', skiprows=cmts, engine='python') 
    n = len(df.columns)
    if n<2:
        print("Need at least two columns in your datafile.  Only found {}, :{}".format(n, df.columns))
        sys.exit(-1)
    if not variables:
        xvar = df.columns[0]
        variables = set(df.columns)
    else:
        variables.intersection(df.columns)
    frames.append(df)

# determine the variable label map
xvarmap = {}

if args.xvar:
    n,l = args.xvar
    xvar = n
    xvarmap[n] = l
else:
    xvarmap[xvar] = xvar

if not xvar in variables:
    print("Independant variable '{}' does not exist, choices are {}'".format(xvar, variables))
    sys.exit(-1)

if not variables:
    print("No variables in common between all files!")
    sys.exit(-1)

varmap = {}

if not args.variables:
    variables.remove(xvar)
    for v in variables:
       varmap[v] = v
else:
    for v in args.variables:
        n,l = v[0]
        if not n in variables:
            print("Variable '{}' does not exist in all the datasets.".format(n))
            sys.exit(-1)
        varmap[n] = l

# pot
#marker_cycler = cycler('marker',['o','s','^','x','p','D','v'])

for v in varmap:
    fig, ax = plt.subplots()
    for i,df in enumerate(frames):
        ax.scatter(df[xvar], df[v], label="Part "+str(i))
    plt.xlabel(xvarmap[xvar])
    plt.ylabel(varmap[v])
    plt.legend()
    plt.show()

