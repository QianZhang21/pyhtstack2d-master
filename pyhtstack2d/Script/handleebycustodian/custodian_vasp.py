#!/usr/bin/env python
from custodian.custodian import Custodian
from custodian.vasp.validators import VaspFilesValidator, VasprunXMLValidator
from custodian.vasp.handlers import VaspErrorHandler, UnconvergedErrorHandler, \
    NonConvergingErrorHandler, FrozenJobErrorHandler, StdErrHandler, \
    WalltimeHandler
from custodian.vasp.jobs import VaspJob

# Define error handlers and validators
handlers = [
    VaspErrorHandler(), FrozenJobErrorHandler(), StdErrHandler(),
    NonConvergingErrorHandler(), WalltimeHandler(),
    UnconvergedErrorHandler()
]
validators = [VaspFilesValidator(), VasprunXMLValidator()]


vasp_cmd = "vasp_std"  # The module name of the VASP executable
max_errors = 3  # The number of errors to fix before giving up

jobs = VaspJob.full_opt_run(vasp_cmd=['mpirun', vasp_cmd], output_file="runlog")
c = Custodian(handlers, jobs, validators=validators, max_errors=max_errors)
c.run()

