#!/bin/bash

project.py --xml beam_on_numu.xml --stage beam_on_numu --clean
project.py --xml beam_on_numu.xml --stage beam_on_numu --submit
project.py --xml beam_off_numu.xml --stage beam_off_numu --clean
project.py --xml beam_off_numu.xml --stage beam_off_numu --submit
project.py --xml bnb_nu_numu.xml --stage bnb_nu_numu --clean
project.py --xml bnb_nu_numu.xml --stage bnb_nu_numu --submit
project.py --xml bnb_nue_numu.xml --stage bnb_nue_numu --clean
project.py --xml bnb_nue_numu.xml --stage bnb_nue_numu --submit
project.py --xml bnb_dirt_numu.xml --stage bnb_dirt_numu --clean
project.py --xml bnb_dirt_numu.xml --stage bnb_dirt_numu --submit
