# Clawpack environment settings
export CLAW=${PROTEUS}"/externalPackages/clawpack"
export CLAWUTIL=${PROTEUS}"/externalPackages/clawpack/clawutil"
export PYCLAW=${PROTEUS}"/externalPackages/clawpack/pyclaw"
export RIEMANN=${PROTEUS}"/externalPackages/clawpack/riemann"
export VISCLAW=${PROTEUS}"/externalPackages/clawpack/visclaw"
export PYTHONPATH="${PROTEUS}/externalPackages/clawpack/visclaw/src/python:${PROTEUS}/externalPackages/clawpack/riemann/src/python:${PROTEUS}/externalPackages/clawpack/pyclaw/src:${PROTEUS}/externalPackages/clawpack/clawutil/src/python:${PYTHONPATH}"
export MATLABPATH=${PROTEUS}"/externalPackages/clawpack/visclaw/src/matlab:${MATLABPATH}"
