import os
import sys
import configparser

def main():
    
    slopes=[0.8]
    imf_type='bi'
    for i in slopes:
        os.system('python newbase/new_base.py {} {}'.format(i, imf_type))
        os.system('conda run -n ambiente2 python2 run_ppxf.py')
        os.system('python determinemass/determine_mass.py {}'.format(i))
        
if __name__ == '__main__':
    main()
