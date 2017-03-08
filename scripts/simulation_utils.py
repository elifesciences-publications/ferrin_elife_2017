import os
import re
import numpy as np
import pandas as pd
from Bio.SeqUtils.CodonUsage import SynonymousCodons
from Bio.SeqUtils import seq1

codonnum = 0
codonDict = dict()
for aa in sorted(SynonymousCodons, key=lambda aa3: seq1(aa3)):
    if aa == 'STOP':
        continue
    for codon in sorted(SynonymousCodons[aa]):
        # these two codons are numbered out of order for consistent notation
        # with Subramaniam et al. Cell 2014
        if codon in ['AGC']:
            codonDict['AGC'] = 59
        elif codon in ['AGT']:
            codonDict['AGT'] = 60
        else:
            codonDict[codon] = codonnum
            codonnum += 1


def get_model(row):
    if row['selpreterm']:
        return 'SAT'
    elif row['5primepreterm']:
        return 'CSAT'
    else:
        return 'TJ'


def get_model_string(row):
    if row['selpreterm']:
        return 'selpreterm'
    elif row['5primepreterm']:
        return '5primepreterm'
    else:
        return 'trafficjam'


def get_parameter(parameter, string):
    try:
        return np.float(
            re.search(parameter + '_([\d\.]+)[_/]', string).groups()[0])
    except AttributeError:
        return 0


def get_nth_line(filehandle, n):
    for linenumber in range(n):
        line = filehandle.readline()
    return line


def get_state(filehandle, state):
    for linenumber in range(n):
        line = filehandle.readline()
    return line


def get_simulation_time(File):
    simulation_time = -1
    filehandle = open(
        File.replace('gene_totetimes.out', 'simulation_parameters.c'))
    for line in filehandle:
        match = re.search('totalSimulatedTime = ([e\d\+\.]+);', line)
        if match:
            simulation_time = np.float(match.groups()[0])
    return simulation_time


def get_stall_strength(row):
    stallstrengthfile1 = open(
        '../annotations/simulations/run2/{:s}_stallstrengthindex_{:d}.tsv'.
        format(get_model_string(row), int(row['stallstrengthindex']))).read()
    stallstrengthfile2 = open(
        '../annotations/simulations/run13/{:s}_stallstrengthindex_{:d}.tsv'.
        format(get_model_string(row), int(row['stallstrengthindex']))).read()
    codon = row['mutant'][:3].upper()
    codonnum = codonDict[codon]
    stallstrength = re.search(
        '\n{:d}\t\w+\t([\w\.]+)\n'.format(codonnum),
        stallstrengthfile1 + stallstrengthfile2).groups()[0]
    return float(stallstrength)


def get_simulation_data(
        runnumber=0,
        ribodensity=False,
        stateoccupancy=False,
        collisionfrequency=False, ):

    rundir = (os.path.abspath(os.path.join(os.getcwd(), os.pardir)) +
              '/rawdata/simulations/run' + str(runnumber) + '/')
    resultdirs = os.listdir(rundir)
    resultdirs = filter(
        lambda x: 'gene_totetimes.out' in os.listdir(rundir + x), resultdirs)
    allFiles = [
        rundir + folder + '/gene_totetimes.out' for folder in resultdirs
    ]

    # find number of proteins produced from YFP0 (1064) and YFP15 (1065) mRNAs.
    data = dict()

    data['proteins_yfp0'] = map(
        lambda x: int(re.search('\n0\W([\d]+)', open(x).read()).groups()[0]),
        allFiles)
    data['proteins_mutant'] = map(
        lambda x: int(re.search('\n1\W([\d]+)', open(x).read()).groups()[0]),
        allFiles)

    # data['simulation_time'] = map( lambda File:
    #                           get_simulation_time(File), allFiles)
    data['files'] = allFiles

    if ribodensity:
        riboFiles = [
            x.replace('gene_totetimes.out', 'gene_ribo_density.out')
            for x in allFiles
        ]

        collisionFiles = [
            x.replace('gene_totetimes.out', 'gene_number_collisions.out')
            for x in allFiles
        ]

        data['ribodensity_yfp0'] = map(
            lambda File: np.array(get_nth_line(open(File), 1).split()),
            riboFiles)
        data['ribodensity_mutant'] = map(
            lambda File: np.array(get_nth_line(open(File), 2).split()),
            riboFiles)

    if collisionfrequency:
        data['nhit5_yfp0'] = map( lambda x: int(
                re.search('\n0\W([-\d]+)\W([-\d]+)', open(x).read()).groups()[0]
            ), collisionFiles)
        data['nhit5_mutant'] = map( lambda x: int(
                re.search('\n1\W([-\d]+)\W([-\d]+)', open(x).read()).groups()[0]
            ), collisionFiles)
        data['nhit3_yfp0'] = map( lambda x: int(
                re.search('\n0\W([-\d]+)\W([-\d]+)', open(x).read()).groups()[1]
            ), collisionFiles)
        data['nhit3_mutant'] = map( lambda x: int(
                re.search('\n1\W([-\d]+)\W([-\d]+)', open(x).read()).groups()[1]
            ), collisionFiles)

    if stateoccupancy:
        stateFiles = [
            x.replace('gene_totetimes.out', 'avg_ribo_tRNA.out')
            for x in allFiles
        ]

        statestrings = {
            'n_rb_ae': 'asiteempty_ribosomes',
            'n_rb_ao': 'asiteoccupied_ribosomes',
            'n_rb_hit3': 'hit3_ribosomes',
            'n_rb_hit5': 'hit5_ribosomes',
            'n_rf': 'free_ribosomes',
        }

        for state, string in statestrings.items():
            data[state] =  map( lambda x: float(
                                                re.search('\n' + string + '\W([\d\.]+)', open(x).read()).groups()[0]
                                            ), stateFiles)

    data['simulation_time'] = map(lambda File: get_simulation_time(File),
                                  allFiles)

    data['mutant'] = map(
        lambda File: re.search('_yfp_(.*)_initiationrate', File).groups()[0],
        allFiles)

    parameterlist = [
        'initiationrate', 'stallstrengthindex', '5primepreterm', 'selpreterm',
        'bkgdpreterm'
    ]

    for parameter in parameterlist:
        data[parameter] = map(lambda File: get_parameter(parameter, File),
                              allFiles)

    data = pd.DataFrame.from_dict(data)
    data['ps_ratio'] = data['proteins_mutant'] / data['proteins_yfp0']
    if runnumber in [2, 13]:
        data['stallstrength'] = data.apply(get_stall_strength, axis=1)
    return data
