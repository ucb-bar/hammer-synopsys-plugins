#!/bin/bash
# Shell script that calls mypy a number of times and then suppresses the output
# if empty so that mypy_with_exclusions.sh can report success.
# Exclude bits known to be problematic to mypy.

set -e

err=0

call_mypy () {
    >&2 echo "Running mypy $*"
    output=$(mypy "$@" | grep -v "python-jsonschema-objects" | grep -v TechJSON | grep -v installs | grep -v tarballs | grep -v ccs | grep -v nldm | grep -v supplies | grep -v lef_file | grep -v qrc_techfile | grep -v serialize | grep -v "hammer_tech.Library" | grep -v milkyway_ | grep -v tluplus | grep -v jsonschema | grep -v "pyyaml/" | grep -v provides | grep -v \"libraries\" || true)
    if [[ ! -z "${output}" ]]; then
        echo "${output}"
        err=1
    fi
}

# Plugins
call_mypy ../../synthesis/dc/__init__.py
call_mypy ../../par/icc/__init__.py
call_mypy ../../drc/icv/__init__.py
call_mypy ../../lvs/icv/__init__.py

exit $err
