import os
import sys
import configparser

def main():
    
    slopes=[1.3, 1.5, 1.8, 2.3, 2.8, 3.3]
    imf_type='bi'
    for i in slopes:
        os.system('python newbase/new_base.py {} {}'.format(i, imf_type))
        os.system('conda run -n ambiente2 python2 run_ppxf.py')
        os.system('python determinemass/determine_mass.py {}'.format(i))
        
if __name__ == '__main__':
    main()
