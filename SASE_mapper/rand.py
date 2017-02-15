from random import random
import numpy as np


def multi_random(range_max, size):
    """
    requires
        range_max: upper bound of range (0 to range_max)
        size: tuple of dimensions
    returns: array of random integers without replacement
    """
	
    sample_arr = np.zeros(size, dtype=int)

    sample_rand = np.random.randint(low=0,high=range_max,size=size)

    for i in xrange(len(sample_arr)):
    	sample = set(sample_rand[i])

    	while len(sample) != size[1]:
    		sample.add(int(random()*range_max))

    	sample_arr[i] = np.array(list(sample))

    return sample_arr


def single_random(range_max, size):
    """
    requires
        range_max: upper bound of range (0 to range_max)
        size: tuple of dimensions
    returns: array of random integers without replacement
    """
	
    sample = set(np.random.randint(low=0,high=range_max,size=size))
    while len(sample) != size:
    	sample.add(int(random()*range_max))

    return np.array(list(sample), dtype=int)


def random_int(range_max, size):
    '''
    Generates random integers without replacement
    
    requires
        range_max: upper bound of range (0 to range_max)
        size: tuple of dimensions
    returns: array of random integers without replacement
    '''
    
    try:
        size[1]
        return multi_random(range_max, size)

    except TypeError:

        return single_random(range_max, size)