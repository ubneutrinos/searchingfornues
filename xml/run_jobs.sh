#!/bin/bash

project.py --xml beam_on.xml --stage beam_on --clean
project.py --xml beam_on.xml --stage beam_on --submit
project.py --xml beam_off.xml --stage beam_off --clean
project.py --xml beam_off.xml --stage beam_off --submit
project.py --xml bnb_nu.xml --stage bnb_nu --clean
project.py --xml bnb_nu.xml --stage bnb_nu --submit
project.py --xml bnb_nue.xml --stage bnb_nue --clean
project.py --xml bnb_nue.xml --stage bnb_nue --submit
project.py --xml bnb_dirt.xml --stage bnb_dirt --clean
project.py --xml bnb_dirt.xml --stage bnb_dirt --submit
