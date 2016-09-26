#pylint: disable=too-many-public-methods, invalid-name, no-self-use
from __future__ import print_function, absolute_import, division

import os
import shutil
import subprocess
import sys
import unittest

from testfixtures.tempdirectory import TempDirectory

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class SnakeFileTest(unittest.TestCase):
    def test_concat_reads(self):
        with TempDirectory() as tmp_dir:
            os.chdir(tmp_dir.path)
            pass