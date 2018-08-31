#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
from os.path import join, isdir
import shutil
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from testfixtures import tempdir

from scripts import unlocker

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

class MockPrompt(object):
    def __init__(self, prompt_return=False):
        self.prompt_was_called = False
        self.prompt_return = prompt_return

    def prompt(self):
        self.prompt_was_called = True
        return self.prompt_return

class UnlockerTest(unittest.TestCase):

    def test_parse_command_line_args(self):
        args = unlocker._parse_command_line_args([])
        self.assertEquals(False, args.dryrun)
        self.assertEquals(False, args.dag)

        args = unlocker._parse_command_line_args(['--dryrun'])
        self.assertEquals(True, args.dryrun)
        self.assertEquals(False, args.dag)

        args = unlocker._parse_command_line_args(['--dag'])
        self.assertEquals(False, args.dryrun)
        self.assertEquals(True, args.dag)

    @tempdir()
    def test_is_dir_locked(self, temp_dir):
        working_dir = temp_dir.path
        self.assertEqual(False, unlocker._is_dir_locked(working_dir))

        snakemake_dir = join(working_dir, '.snakemake')
        os.mkdir(snakemake_dir)
        self.assertEqual(False, unlocker._is_dir_locked(working_dir))

        locks_dir = join(snakemake_dir, 'locks')
        os.mkdir(locks_dir)
        self.assertEqual(False, unlocker._is_dir_locked(working_dir))

        lockfile_path = join(locks_dir, 'foo')
        with open(lockfile_path, 'w') as lockfile:
            print("foo", file=lockfile)
        self.assertEqual(True, unlocker._is_dir_locked(working_dir))

    @tempdir()
    def test_main_skipsUnlockOnDryrun(self, temp_dir):
        working_dir = temp_dir.path
        snakemake_dir = join(working_dir, '.snakemake')
        os.mkdir(snakemake_dir)
        locks_dir = join(snakemake_dir, 'locks')
        os.mkdir(locks_dir)
        lockfile_path = join(locks_dir, 'foo')
        with open(lockfile_path, 'w') as lockfile:
            print("foo", file=lockfile)
        log = StringIO()
        mock_prompt = MockPrompt(True)

        unlocker.main(sys_argv=['--dryrun'],
                      working_dir=working_dir,
                      log=log,
                      prompt_to_unlock=mock_prompt.prompt)

        log_lines = iter(log.getvalue().strip('\n').split('\n'))
        self.assertEqual('skipping lock check (dryrun/dag)', next(log_lines))
        self.assertRegexpMatches(next(log_lines), '===')
        self.assertRaises(StopIteration, next, log_lines)
        self.assertEqual(False, mock_prompt.prompt_was_called)
        self.assertEqual(True, isdir(locks_dir))

    @tempdir()
    def test_main_skipsUnlockOnDag(self, temp_dir):
        working_dir = temp_dir.path
        snakemake_dir = join(working_dir, '.snakemake')
        os.mkdir(snakemake_dir)
        locks_dir = join(snakemake_dir, 'locks')
        os.mkdir(locks_dir)
        lockfile_path = join(locks_dir, 'foo')
        with open(lockfile_path, 'w') as lockfile:
            print("foo", file=lockfile)
        log = StringIO()
        mock_prompt = MockPrompt(True)

        unlocker.main(sys_argv=['--dag'],
                      working_dir=working_dir,
                      log=log,
                      prompt_to_unlock=mock_prompt.prompt)

        log_lines = iter(log.getvalue().strip('\n').split('\n'))
        self.assertEqual('skipping lock check (dryrun/dag)', next(log_lines))
        self.assertRegexpMatches(next(log_lines), '===')
        self.assertRaises(StopIteration, next, log_lines)
        self.assertEqual(False, mock_prompt.prompt_was_called)
        self.assertEqual(True, isdir(locks_dir))

    @tempdir()
    def test_main_skipsUnlockOnCleanDir(self, temp_dir):
        working_dir = temp_dir.path
        locks_dir = join(working_dir, '.snakemake', 'locks')
        _mkdir(locks_dir)
        log = StringIO()
        mock_prompt = MockPrompt(True)

        unlocker.main(sys_argv=[],
                      working_dir=working_dir,
                      log=log,
                      prompt_to_unlock=mock_prompt.prompt)

        log_lines = iter(log.getvalue().strip('\n').split('\n'))
        self.assertEqual('checking for locks: OK', next(log_lines))
        self.assertRegexpMatches(next(log_lines), '===')
        self.assertRaises(StopIteration, next, log_lines)
        self.assertEqual(False, mock_prompt.prompt_was_called)
        self.assertEqual(True, isdir(locks_dir))


    @tempdir()
    def test_main_positivePromptDeletesLockDir(self, temp_dir):
        working_dir = temp_dir.path
        locks_dir = join(working_dir, '.snakemake', 'locks')
        _mkdir(locks_dir)
        lockfile_path = join(locks_dir, 'lockfile')
        with open(lockfile_path, 'w') as lockfile:
            print("foo", file=lockfile)
        log = StringIO()
        mock_prompt = MockPrompt(True)

        exit_code = unlocker.main(sys_argv=[],
                                  working_dir=working_dir,
                                  log=log,
                                  prompt_to_unlock=mock_prompt.prompt)
        log_lines = iter(log.getvalue().strip().split('\n'))
        self.assertRegex(next(log_lines), r'checking .*WARNING: .* locked.')
        self.assertRegex(next(log_lines), r'If you are sure.*')
        self.assertRegex(next(log_lines), r'.* remaining lock.* power loss')
        self.assertRegex(next(log_lines), r'Attempting.*OK')
        self.assertRegexpMatches(next(log_lines), '===')
        self.assertRaises(StopIteration, next, log_lines)
        self.assertEqual(True, mock_prompt.prompt_was_called)
        self.assertEqual(False, isdir(locks_dir))
        self.assertEqual(0, exit_code)

    @tempdir()
    def test_main_negayivePromptStops(self, temp_dir):
        working_dir = temp_dir.path
        locks_dir = join(working_dir, '.snakemake', 'locks')
        _mkdir(locks_dir)
        lockfile_path = join(locks_dir, 'lockfile')
        with open(lockfile_path, 'w') as lockfile:
            print("foo", file=lockfile)
        log = StringIO()
        mock_prompt = MockPrompt(False)

        exit_code = unlocker.main(sys_argv=[],
                                  working_dir=working_dir,
                                  log=log,
                                  prompt_to_unlock=mock_prompt.prompt)
        log_lines = iter(log.getvalue().strip().split('\n'))
        self.assertRegex(next(log_lines), r'checking .*WARNING: .* locked.')
        self.assertRegex(next(log_lines), r'If you are sure.*')
        self.assertRegex(next(log_lines), r'.* remaining lock.* power loss')
        self.assertRegex(next(log_lines), r'Watermelon stopped')
        self.assertRegexpMatches(next(log_lines), '===')
        self.assertRaises(StopIteration, next, log_lines)
        self.assertEqual(True, mock_prompt.prompt_was_called)
        self.assertEqual(True, isdir(locks_dir))
