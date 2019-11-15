#!/bin/bash

project.py --xml beam_off_calo_acpt.xml --stage beam_off_acpt --clean
project.py --xml beam_off_calo_acpt.xml --stage beam_off_acpt --submit
project.py --xml overlay_calo_acpt.xml --stage overlay_acpt --clean
project.py --xml overlay_calo_acpt.xml --stage overlay_acpt --submit
