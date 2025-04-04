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
    def make_fasta(self, input_string):
        with self.temp_fasta.open("w") as ofstream:
            ofstream.write(input_string)
        return Fasta(str(self.temp_fasta), build_index=True)
