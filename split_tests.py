# %load split_tests.py
#!/usr/bin/env python
from typing import List
import re
import glob

filenames: List[str] = glob.glob('tests/testthat_backup/test-*.R')
test_fns: List[str] = []
for fname in filenames:
    s: str = open(fname).read()
    matches: List[str] = re.findall('(test_that\(.*?\\{.*?\}\))', s, flags=re.DOTALL)
    test_fns.extend(matches)

for test_fn in test_fns:
    print(test_fn.split("\n")[0])
    name_matches = re.match('.*["\'](.*)["\'].*', test_fn)
    if name_matches is None: continue
    name = name_matches.group(1)
    name = name.replace(' ', '_')
    print(name)
    open(f'tests/testthat/test-{name}.R', 'w').write(test_fn)
    # TODO: Need to include
    # library(parallel); cores <- detectCores()
    # library(tictoc)
    # Then
    # library(testthat)
    # test_file(f'tests/testthat/test-{name}.R')
