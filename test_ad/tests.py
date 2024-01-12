import pytest
import os
import subprocess


def test_repo_ant40_b2_ss25ka_H():
	subprocess.run (['cp', '../headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -jsf inputs.json', shell = True, check = True)

def test_repo_ant64_b2_future09_ctrl_H():
	subprocess.run (['cp', '../headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_b2_future09_ctrl -dom ant -delta 5.e-2 -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_H():	
	subprocess.run (['cp', '../headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -jsf inputs.json', shell = True, check = True)

def test_repo_ant40_b2_ss25ka_q_geo():
	subprocess.run (['cp', '../headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -iv q_geo -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_q_geo():
        subprocess.run (['cp', '../headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
        subprocess.run('python3 tapenade_config.py -iv q_geo -jsf inputs.json', shell = True, check = True)
