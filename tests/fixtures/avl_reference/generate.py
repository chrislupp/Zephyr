# Project Zephyr
#
# Copyright 2024-2026, Christopher A. Lupp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This work has been cleared for public release, distribution unlimited, case
# number: AFRL-2026-2198.
#
# The views expressed are those of the authors and do not reflect the
# official guidance or position of the United States Government, the
# Department of Defense or of the United States Air Force.
#
# Statement from DoD: The Appearance of external hyperlinks does not
# constitute endorsement by the United States Department of Defense (DoD) of
# the linked websites, of the information, products, or services contained
# therein. The DoD does not exercise any editorial, security, or other
# control over the information you may find at these locations.

"""Regenerate independent references with OptVL; not needed to run Rust tests.

Run with Python containing optvl, from any directory. The wing/tail share one
COMPONENT to match the Rust programmatic fixture's common component index.
AVLs default separate-component core is much wider and changes interference.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import optvl
from optvl import OVLSolver

root = Path(__file__).resolve().parent
library = next(Path(optvl.__file__).parent.glob('libavl*.so'))
parser = argparse.ArgumentParser(__doc__)
parser.add_argument('--checkout', type=Path, required=True, help='source checkout used to build installed OptVL')
checkout = parser.parse_args().checkout
try:
    commit = subprocess.check_output(['git', '-C', str(checkout), 'rev-parse', 'HEAD'], text=True).strip()
except subprocess.CalledProcessError:
    commit = None
result = {'provenance': {'version': optvl.__version__, 'commit': commit,
                        'library_sha256': hashlib.sha256(library.read_bytes()).hexdigest()}, 'cases': []}
for name, alpha, beta in [('rectangular', 4, 2), ('rectangular_claf', 4, 2), ('rectangular_sine', 4, 2), ('wing_tail', 5, 0), ('ucrm_5_8', 2, 0), ('ucrm_17_40', 2, 0)]:
    p = OVLSolver(geo_file=str(root / f'{name}.avl'))
    p.set_variable('alpha', alpha)
    p.set_variable('beta', beta)
    p.execute_run()
    header, forces, derivatives = p.get_header_params(), p.get_total_forces(), p.get_stab_derivs()
    n = p.get_mesh_size()
    case = dict(name=name, alpha_deg=alpha, beta_deg=beta, mach=float(header['mach']),
                iysym=int(header['iysym']), vortices=n,
                coefficients={key: float(forces[key]) for key in ['CL', 'CD', 'Cm']},
                alpha_derivatives={key: float(derivatives[f'd{key}/dalpha']) for key in ['CL', 'CD', 'Cm']},
                circulation=p.get_avl_fort_arr('VRTX_R', 'GAM')[:n].tolist())
    if n < 50:
        case['geometry'] = {key: p.get_avl_fort_arr('VRTX_R', key)[:n, :].tolist()
                            for key in ['RC', 'RV', 'RV1', 'RV2', 'ENC']}
    result['cases'].append(case)
    print(name, case['coefficients'], flush=True)
(root / 'optvl.json').write_text(json.dumps(result, indent=2)+'\n')
