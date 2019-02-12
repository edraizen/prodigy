#!/usr/bin/env python
#
# This code is part of the binding affinity prediction tools distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.
#

"""
Functions to execute freesasa and parse its output.
"""

from __future__ import print_function, division

import os
import subprocess
import sys
import tempfile
from itertools import groupby

try:
    import freesasa
except ImportError as e:
    print('[!] The binding affinity prediction tools require freesasa-python', file=sys.stderr)
    raise ImportError(e)

try:
    from Bio.PDB import PDBParser
    from Bio.PDB import PDBIO, Select
except ImportError as e:
    print('[!] The binding affinity prediction tools require Biopython', file=sys.stderr)
    raise ImportError(e)

from prodigy.config import FREESASA_BIN, FREESASA_PAR
from prodigy.lib.aa_properties import rel_asa

def execute_freesasa(structure, selection=None):
    """
    Runs the freesasa executable on a PDB file.

    You can get the executable from:
        https://github.com/mittinatten/freesasa

    The binding affinity models are calibrated with the parameter
    set for vdW radii used in NACCESS:
        http://www.ncbi.nlm.nih.gov/pubmed/994183
    """
    io = PDBIO()

    param_f = FREESASA_PAR
    # try to get freesasa paths from environment if not use the ones defined in config file

    if 'FREESASA_PAR' in os.environ:
        if os.path.isfile(os.environ['FREESASA_PAR']):
            param_f = os.environ['FREESASA_PAR']

    if param_f is None or not os.path.isfile(param_f):
        param_f = os.path.join(os.path.dirname(os.path.dirname( \
            os.path.abspath(__file__))), "data", "naccess.config")

    # Rewrite PDB using Biopython to have a proper format
    # freesasa is very picky with line width (80 characters or fails!)
    # Select chains if necessary
    class ChainSelector(Select):
        def accept_chain(self, chain):
            if selection and chain.id in selection:
                return 1
            elif not selection:
                return 1
            else:
                return 0

    #Save chains in selection to new PDB file
    #Will remove when freesasa-python has the ability select chain groups
    _pdbf = tempfile.NamedTemporaryFile()
    io.set_structure(structure)
    io.save(_pdbf.name, ChainSelector())

    try:
        # Run freesasa
        classifier = freesasa.Classifier(param_f)
        structure = freesasa.Structure(_pdbf.name, classifier)
        result = freesasa.calc(structure)
    except (KeyboardInterrupt, SystemExit):
        #freesasa will throw a regular exception
        raise
    except (AssertionError, IOError, Exception) as e:
        print('[!] freesasa did not run successfully', file=sys.stderr)
        print(e, file=sys.stderr)
        raise Exception(e)

    # Parse freesasa results
    asa, rsa = parse_freesasa_output(result, structure)

    _pdbf.close()

    return asa, rsa

def parse_freesasa_output(result, structure):
    """
    Returns per-residue relative accessibility of side-chain and main-chain
    atoms as calculated by freesasa.
    """

    asa_data, rsa_data = {}, {}

    _rsa = rel_asa
    _bb = set(('CA', 'C', 'N', 'O'))

    for res, atoms in groupby(range(result.nAtoms()), key=structure.residueNumber):
        res_id = None
        asa_mc, asa_sc, total_asa = 0, 0, 0
        for i in atoms:
            at_id = (structure.chainLabel(i), structure.residueName(i),
                structure.residueNumber(i), structure.atomName(i))
            if res_id is None:
                res_id = at_id[:-1]
            asa = result.atomArea(i)
            # if structure.atomName(i) in _bb:
            #     asa_mc += asa
            # else:
            #     asa_sc += asa
            total_asa += asa
            asa_data[at_id] = asa

        rsa_data[res_id] = total_asa / _rsa['total'][structure.residueName(i)]

    return asa_data, rsa_data
