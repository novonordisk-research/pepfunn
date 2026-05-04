"""
PepFunNN: Protocols for the analysis of peptides using cheminformatics and bioinformatics tools
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__email__ = "raoc@novonordisk.com"

########################################################################################
# Modules to import
########################################################################################

import os
import time
import unittest

import pepfunn.similarity as _sim
from pepfunn.similarity import pepDescriptors, monomerFP
from pepfunn.sequence import SequenceConstants

##########################################################################
# Helpers
##########################################################################

_SEQS = [
    'AFTGYW', 'AGTGYL', 'LLSHYTSY', 'NPVVHFFKNIVTPRTPPPSQ',
    'W-W-S-E-V-N-R-A-E-F', 'K-T-E-E-I-S-E-V-N-I-V-A-E-F',
    'K-Aib-M-P', 'S-A-Aib-P', 'A-F-T-G-Y-W', 'L-L-S-H-Y-T-S-Y',
]

def _prop_file():
    module_dir = os.path.dirname(os.path.abspath(_sim.__file__))
    return os.path.join(module_dir, SequenceConstants.def_path, SequenceConstants.def_property)

########################################################################################

class TestMonomerFPPerformance(unittest.TestCase):

    def test_prop_lookup_speedup(self):
        """
        Property lookup: one file-read + DataFrame scan per call vs module-level dict cache.
        Asserts the cached path is at least 50x faster.
        """
        file_path = _prop_file()
        monomers = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        prop_list = ['heavy', 'nrot', 'hacc', 'hdon', 'nhet', 'tpsa', 'mw']
        n = 200

        t0 = time.perf_counter()
        for _ in range(n):
            with open(file_path, 'r') as fh:
                df = pepDescriptors.get_properties(fh)
            for mon in monomers:
                for prop in prop_list:
                    df.loc[df['name'] == mon, prop].item()
        slow = time.perf_counter() - t0

        _sim._prop_cache.clear()
        monomerFP('AFTGYW')
        lookup = _sim._prop_cache[file_path]

        t0 = time.perf_counter()
        for _ in range(n):
            for mon in monomers:
                for prop in prop_list:
                    lookup[mon][prop]
        fast = time.perf_counter() - t0

        speedup = slow / fast
        print(f"\n  naive {n}x:  {slow:.3f}s")
        print(f"  cached {n}x: {fast:.4f}s")
        print(f"  speedup:     {speedup:.0f}x")

        self.assertGreater(speedup, 50)

    ##########################################################################

    def test_batch_monomerFP(self):
        """
        End-to-end timing for monomerFP over all pairs in _SEQS.
        Printed for use as benchmark evidence; asserts reasonable wall time.
        """
        pairs = [(s1, s2) for i, s1 in enumerate(_SEQS) for s2 in _SEQS[i:]]
        _sim._prop_cache.clear()

        t0 = time.perf_counter()
        for s1, s2 in pairs:
            monomerFP(s1)
            monomerFP(s2)
        elapsed = time.perf_counter() - t0

        ms_per_pair = 1000 * elapsed / len(pairs)
        print(f"\n  {len(pairs)} pairs: {elapsed:.3f}s ({ms_per_pair:.2f}ms/pair)")

        self.assertLess(elapsed, 120)

########################################################################################

if __name__ == "__main__":
    unittest.main()
