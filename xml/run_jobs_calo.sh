#!/bin/bash

project.py --xml beam_on_calo.xml --stage beam_on --clean
project.py --xml beam_on_calo.xml --stage beam_on --submit
project.py --xml beam_off_calo.xml --stage beam_off --clean
project.py --xml beam_off_calo.xml --stage beam_off --submit
project.py --xml bnb_nu_calo.xml --stage bnb_nu --clean
project.py --xml bnb_nu_calo.xml --stage bnb_nu --submit
project.py --xml bnb_dirt_calo.xml --stage bnb_dirt --clean
project.py --xml bnb_dirt_calo.xml --stage bnb_dirt --submit
