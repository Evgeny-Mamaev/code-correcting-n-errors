import unittest

from linearcode import ErrorCorrectingCode
from fileutils import remove_files


class BuildCodeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        remove_files()

    def test_capability(self):
        ErrorCorrectingCode(r=19, n=23, t=4, channel_error_probability=0.01)
