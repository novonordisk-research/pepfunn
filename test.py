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

from unittest import TestLoader, TestResult
from pathlib import Path

##########################################################################
# Function
##########################################################################

def run_tests():
    """
    Function to run all the available unittests
    """

    test_loader = TestLoader()
    test_result = TestResult()

    test_directory = str(Path(__file__).resolve().parent / 'tests')

    test_suite = test_loader.discover(test_directory, pattern='test_*.py')
    test_suite.run(result=test_result)

    if test_result.wasSuccessful():
        print("\n################################\nAll of the tests were successful!\n################################")
        exit(0)
    else:
        print("\nFailures:")
        for i,fail in enumerate(test_result.failures):
            print(f"{i+1}. {fail[0]}\n{fail[1]}")
        if test_result.errors:
            print("\nErrors:")
            for i,error in enumerate(test_result.errors):
                print(error)
        exit(-1)

if __name__ == '__main__':
    run_tests()