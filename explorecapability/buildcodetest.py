import unittest

from linearcode import ErrorCorrectingCode
from test.testutils import remove_files


class BuildCodeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        remove_files()

    def test_capability(self):
        ErrorCorrectingCode(r=23, n=28, t=5, channel_error_probability=0.01)
