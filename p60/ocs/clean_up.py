#! /usr/local/bin/python

import ocs, sys 

datatran = ocs.ocs_ooriented()
datatran.clean_up(sys.argv[1], number=int(sys.argv[2]), wait=int(sys.argv[3])) 


