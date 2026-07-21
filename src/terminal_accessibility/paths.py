"""Input-path expansion. Stdlib only -- no scientific deps, so it can be
imported to gather work without pulling in numpy/biotite."""

import glob
import os


def gather_paths(inputs):
    """Expand files, globs, and directories into a sorted unique list of .pdb/.cif paths."""
    paths = []
    for item in inputs:
        if os.path.isdir(item):
            for root, _, files in os.walk(item):
                paths.extend(os.path.join(root, f) for f in files if f.lower().endswith((".pdb", ".cif")))
        else:
            paths.extend(glob.glob(item))
    return sorted(set(paths))
