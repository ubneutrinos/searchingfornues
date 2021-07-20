#!/bin/bash

#project.py --xml grid_run1_detsys_calo.xml --stage cv --clean
#project.py --xml grid_run1_detsys_calo.xml --stage cv --submit

#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_x --clean
#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_x --submit

#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_yz --clean
#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_yz --submit

#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_theta_xz --clean
#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_theta_xz --submit

#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_theta_yz --clean
#project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_theta_yz --submit

project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_dedx --clean
project.py --xml grid_run1_detsys_calo.xml --stage wire_mod_dedx --submit

project.py --xml grid_run1_detsys_calo.xml --stage sce --clean
project.py --xml grid_run1_detsys_calo.xml --stage sce --submit

project.py --xml grid_run1_detsys_calo.xml --stage recomb2 --clean
project.py --xml grid_run1_detsys_calo.xml --stage recomb2 --submit
