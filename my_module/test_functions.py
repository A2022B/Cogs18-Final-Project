"""Test for my functions.

Note: because these are 'empty' functions (return None), here we just test
  that the functions execute, and return None, as expected.
"""

from functions import astroImage, fourimgshow, basic_stats

def test_astroImage():
    assert callable(astroImage)


def test_fourimgshow():
    assert callable(fourimgshow)

    
def test_basic_stats():
    assert callable(basic_stats) 

#Notes to grader: The problem with my functions are that the outputs are either image or NoneType. I don't see any clear way to test my functions besides "assert callable()". 
