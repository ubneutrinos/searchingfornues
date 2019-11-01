#!/bin/bash

project.py --xml bnb_nu_calo_shr.xml --stage bnb_nu_shr --clean
project.py --xml bnb_nu_calo_shr.xml --stage bnb_nu_shr --submit
project.py --xml bnb_nue_calo_shr.xml --stage bnb_nue_shr --clean
project.py --xml bnb_nue_calo_shr.xml --stage bnb_nue_shr --submit
