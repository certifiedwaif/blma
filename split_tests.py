# %load split_tests.py
#!/usr/bin/env python
from typing import List
import re
import glob
from jinja2 import Template

filenames: List[str] = glob.glob('tests/testthat_backup/test-*.R')
test_fns: List[str] = []
for fname in filenames:
    s: str = open(fname).read()
    matches: List[str] = re.findall('(test_that\(.*?\\{.*?\}\))', s, flags=re.DOTALL)
    test_fns.extend(matches)

preamble: str = """
library(testthat)
library(blma)
library(tictoc)
library(parallel)

cores <- detectCores()

"""
names = []
for test_fn in test_fns:
    # print(test_fn.split("\n")[0])
    name_matches = re.match('.*["\'](.*)["\'].*', test_fn)
    if name_matches is None: continue
    name = name_matches.group(1)
    name = name.replace(' ', '_')
    names.append(name)
    # print(name)
    open(f'tests/testthat/test-{name}.R', 'w').write(preamble + test_fn)
    # TODO: Need to include
    # library(parallel); cores <- detectCores()
    # library(tictoc)
    # Then
    # library(testthat)
    # test_file(f'tests/testthat/test-{name}.R')
job_yaml_template:str = Template("""apiVersion: batch/v1
kind: Job
metadata:
  name: blma
spec:
  template:
    parallelism: {{ len(names) }}
    spec:
      containers: {% for name in names %}
      - name: {{ name.lower().replace('_', '-') }}
        image: certifiedwaif/blma
        command: ["R", "-e", "testthat::test_file('tests/testthat/test-{{ name }}.R')"]
    {% endfor %}
      restartPolicy: Never
  backoffLimit: 4""")
job_yaml: str = job_yaml_template.render(names=names)
open('job.yaml', 'w+').write(job_yaml)
