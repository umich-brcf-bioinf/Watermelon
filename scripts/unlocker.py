#!/usr/bin/env python
from __future__ import print_function, absolute_import, division

from os import listdir, getcwd
from os.path import isdir, isfile, join
import shutil
import sys


_LOCK_WARNING = \
'''WARNING: It looks like the directory is locked.
If you are sure that no other instances of snakemake are running on this directory,
the remaining lock was likely caused by a kill signal or a power loss.
'''

_HEADER_RULE = '=' * 70 + '\n'

def _get_locks_dir(working_dir):
    return join(working_dir, '.snakemake', 'locks')

def _is_dir_locked(working_dir):
    locks_dir = join(working_dir, '.snakemake', 'locks')
    return isdir(locks_dir) and len(listdir(locks_dir)) > 0

def _prompt_to_unlock():
    value = input('Do you want to remove the lock and proceed? (yes/no): ')
    return value.lower().strip() == 'yes'

def _clear_locks(working_dir, log, prompt_to_unlock):
    exit_code = 1
    if prompt_to_unlock():
        log.write('Attempting to clear locks: ')
        shutil.rmtree(_get_locks_dir(working_dir))
        if _is_dir_locked(working_dir):
            raise ValueError('ERROR: Could not clear locks')
        log.write('OK\n')
        exit_code = 0
    else:
        log.write('Watermelon stopped to manually review/correct lock files.\n')
        exit_code = 1
    return exit_code

def main(working_dir=getcwd(),
         log=sys.stderr,
         prompt_to_unlock=_prompt_to_unlock):
    log.write('checking for locks: ')
    exit_code = 0
    if _is_dir_locked(working_dir):
        log.write(_LOCK_WARNING)
        exit_code = _clear_locks(working_dir, log, prompt_to_unlock)
    else:
        log.write('OK\n')
    log.write(_HEADER_RULE)
    return exit_code

if __name__ == '__main__':
    exit_code = main()
    exit(exit_code)