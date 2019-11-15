#!/bin/bash

project.py --xml bnb_ccpi0_calo_shr.xml --stage bnb_ccpi0_shr --clean
project.py --xml bnb_ccpi0_calo_shr.xml --stage bnb_ccpi0_shr --submit
project.py --xml bnb_nue_calo_shr.xml --stage bnb_nue_shr --clean
project.py --xml bnb_nue_calo_shr.xml --stage bnb_nue_shr --submit
project.py --xml bnb_ncpi0_calo_shr.xml --stage bnb_ncpi0_shr --clean
project.py --xml bnb_ncpi0_calo_shr.xml --stage bnb_ncpi0_shr --submit
