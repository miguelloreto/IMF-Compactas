#!/usr/bin/env python
from __future__ import print_function
from ppxf_population_gas_example import ppxf_population_gas_example
if __name__ == '__main__':
    ppxf_population_gas_example('0446-51899-0541.cxt', 'grid_in.dat', 3900, ngas = 2, m_stars = 2, v_stars = 0.0, s_stars = 100, m_gas1 = 2, v_gas1 = 0.0, s_gas1 = 100, v_gas2 = 0.0, s_gas2 = 50.0)
