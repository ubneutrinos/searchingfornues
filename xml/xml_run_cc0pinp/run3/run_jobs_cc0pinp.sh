#!/bin/bash

project.py --xml beam_on_cc0pinp.xml --stage run3_beam_on_cc0pinp --clean
project.py --xml beam_on_cc0pinp.xml --stage run3_beam_on_cc0pinp --submit
project.py --xml beam_off_cc0pinp.xml --stage run3_beam_off_cc0pinp --clean
project.py --xml beam_off_cc0pinp.xml --stage run3_beam_off_cc0pinp --submit
project.py --xml bnb_nu_cc0pinp.xml --stage run3_bnb_nu_cc0pinp --clean
project.py --xml bnb_nu_cc0pinp.xml --stage run3_bnb_nu_cc0pinp --submit
project.py --xml bnb_nue_cc0pinp.xml --stage run3_bnb_nue_cc0pinp --clean
project.py --xml bnb_nue_cc0pinp.xml --stage run3_bnb_nue_cc0pinp --submit
project.py --xml bnb_dirt_cc0pinp.xml --stage run3_bnb_dirt_cc0pinp --clean
project.py --xml bnb_dirt_cc0pinp.xml --stage run3_bnb_dirt_cc0pinp --submit

project.py --xml bnb_elee_low_cc0pinp.xml --stage run3_bnb_elee_low_cc0pinp --clean
project.py --xml bnb_elee_low_cc0pinp.xml --stage run3_bnb_elee_low_cc0pinp --submit
project.py --xml bnb_elee_high_cc0pinp.xml --stage run3_bnb_elee_high_cc0pinp --clean
project.py --xml bnb_elee_high_cc0pinp.xml --stage run3_bnb_elee_high_cc0pinp --submit
