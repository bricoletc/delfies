from dataclasses import dataclass
from pathlib import Path
from random import choice as rand_choice
from tempfile import NamedTemporaryFile

import pytest
from pysam import AlignedSegment, AlignmentFile
from pysam import index as pysam_index
from pysam import qualitystring_to_array

from delfies import BreakpointType, Orientation
from delfies.breakpoint_foci import BreakpointDetectionParams, find_breakpoint_foci
from delfies.interval_utils import Interval
from delfies.SAM_utils import DEFAULT_MIN_MAPQ, DEFAULT_READ_FILTER_FLAG
from delfies.seq_utils import randomly_substitute, rev_comp
from delfies.telomere_utils import TELOMERE_SEQS

DEFAULT_NON_TELO_SEQ = "TAACCC"
DEFAULT_TELO_SEQ = TELOMERE_SEQS["Nematoda"][Orientation.forward]
DEFAULT_CHROM = "chr1"
DEFAULT_CHROM_LENGTH = 10e03
DEFAULT_MIN_TELO_ARRAY_SIZE = 10
DEFAULT_NUM_TELO_CONTAINING_READS = 150
EXPECTED_BREAKPOINT_POSITION = 1000


@dataclass
class TeloContainingReadGenerator:
    read_template: str = DEFAULT_NON_TELO_SEQ * 100 + DEFAULT_TELO_SEQ * 100
    num_reads: int = 200
    num_telo_containing_reads: int = DEFAULT_NUM_TELO_CONTAINING_READS
    telo_start_pos_in_genome: int = EXPECTED_BREAKPOINT_POSITION
    telo_start_pos_in_read: int = len(DEFAULT_NON_TELO_SEQ) * 100
    telo_array_size: int = len(DEFAULT_TELO_SEQ) * DEFAULT_MIN_TELO_ARRAY_SIZE * 2
    num_mutations: int = 0
    BAM_file: str = None

    def _add_telo_containing_softclipped_sequence(self, read: AlignedSegment) -> None:
        telo_start_pos = self.telo_start_pos_in_read
        index_start = rand_choice(range(telo_start_pos))
        non_telo_seq = self.read_template[index_start:telo_start_pos]
        telo_seq = self.read_template[
            telo_start_pos : telo_start_pos + self.telo_array_size
        ]
        read.reference_start = self.telo_start_pos_in_genome - len(non_telo_seq)
        read.cigar = ((0, len(non_telo_seq)), (4, len(telo_seq)))
        read.query_sequence = randomly_substitute(
            non_telo_seq, self.num_mutations
        ) + randomly_substitute(telo_seq, self.num_mutations)
        read.template_length = len(read.query_sequence)
        read.query_qualities = qualitystring_to_array("<" * len(read.query_sequence))

    def _add_nontelo_containing_fully_aligned_sequence(
        self, read: AlignedSegment
    ) -> None:
        telo_start_pos = self.telo_start_pos_in_read
        index_start = rand_choice(range(telo_start_pos))
        non_telo_seq_1 = self.read_template[index_start:telo_start_pos]
        non_telo_seq_2 = self.read_template[
            telo_start_pos - self.telo_array_size : telo_start_pos
        ]
        read.reference_start = self.telo_start_pos_in_genome - len(non_telo_seq_1)
        read.query_sequence = randomly_substitute(
            non_telo_seq_1, self.num_mutations
        ) + randomly_substitute(non_telo_seq_2, self.num_mutations)
        read.cigar = ((0, len(read.query_sequence)),)
        read.template_length = len(read.query_sequence)
        read.query_qualities = qualitystring_to_array("<" * len(read.query_sequence))

    def write_BAM(
        self,
        bitwise_flag: int = 0,
        mapping_quality: int = DEFAULT_MIN_MAPQ,
    ) -> str:
        header = {
            "HD": {"VN": "1.0"},
            "SQ": [{"LN": DEFAULT_CHROM_LENGTH, "SN": DEFAULT_CHROM}],
        }
        BAM_file = NamedTemporaryFile()
        with AlignmentFile(BAM_file.name, "wb", header=header) as outf:
            all_reads = list()
            for i in range(self.num_reads):
                read = AlignedSegment()
                read.reference_id = 0
                read.query_name = f"simulated_read_{i}"
                if i < self.num_telo_containing_reads:
                    self._add_telo_containing_softclipped_sequence(read)
                else:
                    self._add_nontelo_containing_fully_aligned_sequence(read)
                read.flag = bitwise_flag
                read.mapping_quality = mapping_quality
                all_reads.append(read)
            all_reads = sorted(all_reads, key=lambda el: el.reference_start)
            for read in all_reads:
                outf.write(read)
        self.BAM_file = BAM_file
        pysam_index(self.BAM_file.name)
        return self.BAM_file.name

    def __del__(self):
        if self.BAM_file is not None:
            self.BAM_file.close()
            Path(f"{self.BAM_file.name}.bai").unlink()


@pytest.fixture
def read_generator():
    return TeloContainingReadGenerator()


@pytest.fixture
def genome_interval():
    return Interval(DEFAULT_CHROM, 0, DEFAULT_CHROM_LENGTH)


@pytest.fixture
def detection_params():
    telomere_seqs = {
        Orientation.forward: DEFAULT_TELO_SEQ,
        Orientation.reverse: rev_comp(DEFAULT_TELO_SEQ),
    }
    return BreakpointDetectionParams(
        bam_fname="missing",
        telomere_seqs=telomere_seqs,
        telo_array_size=DEFAULT_MIN_TELO_ARRAY_SIZE,
        max_edit_distance=0,
        clustering_threshold=1,
        min_mapq=DEFAULT_MIN_MAPQ,
        read_filter_flag=DEFAULT_READ_FILTER_FLAG,
        min_supporting_reads=10,
        breakpoint_type=BreakpointType.S2G,
    )


def test_forward_breakpoint_S2G(read_generator, genome_interval, detection_params):
    detection_params.bam_fname = read_generator.write_BAM()
    foci = find_breakpoint_foci(detection_params, genome_interval)
    filtered_foci = [
        elem for elem in foci if elem.start == EXPECTED_BREAKPOINT_POSITION
    ]
    assert len(filtered_foci) == 1
    breakpoint_focus = filtered_foci[0]
    assert breakpoint_focus.contig == DEFAULT_CHROM
    assert breakpoint_focus.start == EXPECTED_BREAKPOINT_POSITION
    assert breakpoint_focus.end == EXPECTED_BREAKPOINT_POSITION + 1
    assert breakpoint_focus.breakpoint_type == str(BreakpointType.S2G)
    assert (
        breakpoint_focus.num_supporting_reads__forward
        == DEFAULT_NUM_TELO_CONTAINING_READS
    )
