from datasci import Tents
from edlib import align as edlib_align
from pyfastx import Fasta

from delfies import Orientation
from delfies.SAM_utils import SoftclippedRead, SoftclippedReads
from delfies.seq_utils import rev_comp

TELOMERE_SEQS = {
    "Nematoda": {Orientation.forward: "TTAGGC", Orientation.reverse: "GCCTAA"}
}


def has_softclipped_telo_array(
    read: SoftclippedRead,
    orientation: Orientation,
    telomere_seqs,
    min_telo_array_size: int,
    max_edit_distance: int,
) -> bool:
    """
    Note: we allow for the softclipped telo array to start with any cyclic shift
    of the telomeric repeat unit.
    """
    telo_unit = telomere_seqs[orientation]
    searched_telo_array = telo_unit * min_telo_array_size
    subseq_clip_end = len(searched_telo_array) + len(telo_unit)
    if orientation is Orientation.forward:
        end = read.sc_query + subseq_clip_end
        subseq = read.sequence[read.sc_query : end]
    else:
        start = max(read.sc_query + 1 - subseq_clip_end, 0)
        subseq = read.sequence[start : read.sc_query + 1]
    result = edlib_align(
        searched_telo_array, subseq, mode="HW", task="distance", k=max_edit_distance
    )
    found_telo_array = result["editDistance"] != -1
    return found_telo_array


def setup_telomere_tents() -> Tents:
    tents_header = [
        "contig",
        "start",
        "end",
        "template_nucleotide",
        "added_telomere",
        "sample_name",
    ]
    tents = Tents(header=tents_header)
    return tents


def record_telomere_additions(
    reads_for_telomere_additions: SoftclippedReads,
    telomere_unit_length: int,
    genome_fname: str,
    sample_name: str,
) -> Tents:
    genome = Fasta(genome_fname)
    telomere_tents = setup_telomere_tents()
    for read in reads_for_telomere_additions:
        sc_query = read.sc_query
        if read.orientation is Orientation.forward:
            added_telomere = read.sequence[sc_query : sc_query + telomere_unit_length]
            ref_position_one_based = read.sc_ref
            template_nucleotide = genome.fetch(
                read.reference_name, (ref_position_one_based, ref_position_one_based)
            )
        else:
            added_telomere = rev_comp(
                read.sequence[sc_query + 1 - telomere_unit_length : sc_query + 1]
            )
            ref_position_one_based = read.sc_ref + 2
            template_nucleotide = rev_comp(
                genome.fetch(
                    read.reference_name,
                    (ref_position_one_based, ref_position_one_based),
                )
            )
        new_tent = telomere_tents.new()
        new_tent.update(
            contig=read.reference_name,
            start=ref_position_one_based - 1,
            end=ref_position_one_based,
            template_nucleotide=template_nucleotide,
            added_telomere=added_telomere,
            sample_name=sample_name,
        )
        telomere_tents.add(new_tent)
    return telomere_tents
