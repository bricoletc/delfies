from pathlib import Path
from tempfile import TemporaryDirectory

from pyfastx import Fasta


class ClassWithTempFasta:
    @classmethod
    def setup_class(cls):
        cls.temp_dir = TemporaryDirectory()
        cls.temp_fasta = Path(cls.temp_dir.name) / "temp.fasta"

    @classmethod
    def teardown_class(cls):
        cls.temp_dir.cleanup()

    @classmethod
    def make_fasta(self, input_string, as_handle=False):
        with self.temp_fasta.open("w") as ofstream:
            ofstream.write(input_string)
        if as_handle:
            return Fasta(str(self.temp_fasta), build_index=True)
        else:
            return str(self.temp_fasta)
