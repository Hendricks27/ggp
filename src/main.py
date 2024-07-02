import os
import re
import sys
import gzip
import argparse


class GraphicalFragmentAssemblySegmentForPhasing(object):

    def __init__(self):
        self.clear()

    def clear(self):
        self._segment_length = {}
        self._maternal_walk_list = {}
        self._paternal_walk_list = {}

        self._maternal_unique_segments = {}
        self._paternal_unique_segments = {}

    @staticmethod
    def parse_walk(walk_str):
        # Need to consider orientation
        # Example:
        # >1>2>3>4>5<6>7>8>9>10>11>12>13>14>15>16>17
        res = []

        this_segments = ["", None]
        for c in walk_str:
            if c in "<>":
                if this_segments[0] != "":
                    res.append(this_segments)
                    this_segments = ["", None]

                orientation = 1 if c == ">" else -1
                this_segments[1] = orientation
            else:
                this_segments[0] += c

        if this_segments[0] != "":
            res.append(this_segments)
        assert len(res) == walk_str.count(">") + walk_str.count("<")
        return res

    def get_unique_segments(self):

        # Find unique segments
        psegs = set()
        for k, v in self._paternal_walk_list.items():
            for seg, ori in v:
                psegs.add(seg)

        msegs = set()
        for k, v in self._maternal_walk_list.items():
            for seg, ori in v:
                msegs.add(seg)

        common_segs = psegs.intersection(msegs)

        mat_unique_segments = msegs - common_segs
        pat_unique_segments = psegs - common_segs

        # Key: segment ID
        # Value: [parent_info, child_info]
        # Parent/child info: [segment, edge_orientation]

        mat_phasing_info = {}
        pat_phasing_info = {}

        for seg in mat_unique_segments:
            mat_phasing_info[seg] = [None, None]

        for seg in pat_unique_segments:
            pat_phasing_info[seg] = [None, None]

        # This check only ensure the orientation is consistent for the unique segments.
        # For the common segment, there might be unique edges. (TODO list)
        for paternal_walk_meta, paternal_walk_segments in self._paternal_walk_list.items():
            for i in range(0, len(paternal_walk_segments) - 1):
                this_seg, this_ori = paternal_walk_segments[i]
                next_seg, next_ori = paternal_walk_segments[i + 1]

                # segment orientation is either 1 or -1
                # edge orientation is also either 1 or -1
                # edge orientation 1 means the next segment is in the same direction as the current segment
                # edge orientation -1 means the next segment is in the opposite direction as the current segment
                edge_orientation = this_ori * next_ori

                if this_seg in pat_unique_segments:
                    pat_phasing_info[this_seg][1] = (next_seg, edge_orientation)
                if next_seg in pat_unique_segments:
                    pat_phasing_info[next_seg][0] = (this_seg, edge_orientation)

        for maternal_walk_meta, maternal_walk_segments in self._maternal_walk_list.items():
            for i in range(0, len(maternal_walk_segments) - 1):
                this_seg, this_ori = maternal_walk_segments[i]
                next_seg, next_ori = maternal_walk_segments[i + 1]

                edge_orientation = this_ori * next_ori

                if this_seg in mat_unique_segments:
                    mat_phasing_info[this_seg][1] = (next_seg, edge_orientation)
                if next_seg in mat_unique_segments:
                    mat_phasing_info[next_seg][0] = (this_seg, edge_orientation)

        self._maternal_unique_segments = mat_phasing_info
        self._paternal_unique_segments = pat_phasing_info

        return

    def parse(self, gfa_file, maternal_walk_ids, paternal_walk_ids, ignore_sex_chr=False):
        # maternal_walk_list = {}
        # paternal_walk_list = {}

        never_found_correct_walk = True
        with open(gfa_file) as gfa_fh:
            for l in gfa_fh:
                if l[0] not in "SW":
                    continue

                l = l.strip().split("\t")

                if l[0] == "S":
                    rt, segID, seq, *tags = l
                    # sequence, tag, links, parent_links_count
                    self._segment_length[segID] = len(seq)

                if l[0] == "W":
                    rt, sampleID, parental_index, chrom, start, end, walk_str = l

                    if chrom in ["chrX", "chrY"] and ignore_sex_chr:
                        continue

                    parental_index = int(parental_index)
                    start = int(start)
                    end = int(end)

                    is_maternal_walk = (sampleID, parental_index) == (maternal_walk_ids[0], maternal_walk_ids[1])
                    is_paternal_walk = (sampleID, parental_index) == (paternal_walk_ids[0], paternal_walk_ids[1])

                    if not (is_maternal_walk or is_paternal_walk):
                        continue

                    never_found_correct_walk = False

                    # print(sampleID, parental_index, chrom, start, end, walk_str[:100])

                    walk_list = self.parse_walk(walk_str)
                    key = (chrom, start, end)
                    if is_maternal_walk:
                        self._maternal_walk_list[key] = walk_list
                        # print("Maternal", key, walk_list[:3])
                    if is_paternal_walk:
                        self._paternal_walk_list[key] = walk_list
                        # print("Paternal", key, walk_list[:3])

        if never_found_correct_walk:
            raise RuntimeError(
                f"Finished parsing GFA file, but never found the correct walk. \nMaternal walk provided: {maternal_walk}, \nPaternal walk provided: {paternal_walk}")

        self.get_unique_segments()
        return

    def get_sequence_length_by_segment_ID(self, segment_ID):
        return self._segment_length[segment_ID]


class GraphAlignmentFormat(object):
    cs_tag_regex = re.compile(r"(:\d+|[+-]\w+|\*\w{2})")

    def __init__(self, fp):
        self._fp = fp

    @staticmethod
    def parse_path(path_str):
        # >1>2>3>4>5<6>7>8>9>10>11>12>13>14>15>16>17
        res = []

        this_segments = ["", None]
        for c in path_str:
            if c in "<>":
                if this_segments[0] != "":
                    res.append(this_segments)
                    this_segments = ["", None]

                orientation = 1 if c == ">" else -1
                this_segments[1] = orientation
            else:
                this_segments[0] += c

        if this_segments[0] != "":
            res.append(this_segments)
        assert len(res) == path_str.count(">") + path_str.count("<")
        return res

    def parse_cs_tag(self, cs_tag_str):
        # :30*AG+TTTT-AAAA

        res = []
        reconstructed = ""
        for m in self.cs_tag_regex.findall(cs_tag_str):
            # res.append(m)
            reconstructed += m

            if m[0] == ":":
                res.append(("M", int(m[1:])))
            elif m[0] == "+":
                res.append(("I", len(m[1:])))
            elif m[0] == "-":
                res.append(("D", len(m[1:])))
            elif m[0] == "*":
                res.append(("X", m[1:]))
            else:
                raise RuntimeError
        assert reconstructed == cs_tag_str

        return res

    def get_alignments(self):
        # Assume sorted by query name
        count = 0
        x = 0
        same_read_pair = ["", [], []]

        with open(self._fp) as gaf_fh:
            for l in gaf_fh:
                l = l.strip().split("\t")
                query_name = l[0]

                regular_fields = l[:12]
                tag_list = l[12:]
                tags = {}

                skip = False
                for i in [1, 2, 3, 6, 7, 8, 9, 10, 11]:
                    if regular_fields[i] == "*":
                        skip = True
                        break
                if skip:
                    continue

                for t in tag_list:
                    # print(t)
                    tn, tt, tv = t.split(":", 2)
                    if tt == "i":
                        tags[tn] = int(tv)
                    elif tt == "f":
                        tags[tn] = float(tv)
                    else:
                        tags[tn] = tv

                aln = (regular_fields, tags)
                r1flag = "fn" in tags.keys()
                r2flag = "fp" in tags.keys()

                # print(r1flag, r2flag, tags.keys())
                assert r1flag or r2flag
                # continue

                if len(same_read_pair[0]) == 0 or query_name == same_read_pair[0]:
                    same_read_pair[0] = query_name
                    if r1flag:
                        same_read_pair[1].append(aln)
                    if r2flag:
                        same_read_pair[2].append(aln)
                    continue

                x += 1
                # if x > 100:
                #    break

                # print(same_read_pair)

                # Use "at" tag to describe alignment type
                # First digit: 0 for unique alignment, 1 for multimapping
                # Second digit: 0 for primary alignment, 1 for secondary alignment

                for ri in [1, 2]:
                    read_alignment = same_read_pair[ri]
                    if len(read_alignment) == 0:
                        same_read_pair[ri] = None
                        continue

                    if len(read_alignment) == 1:
                        read_alignment[0][1]["at"] = "00"
                        same_read_pair[ri] = read_alignment[0]
                        continue

                    best_alignment = read_alignment[0]
                    best_score = best_alignment[1]["AS"]
                    for a in read_alignment[1:]:
                        # print(a[1]["AS"], best_score)
                        if a[1]["AS"] > best_score:
                            best_alignment = a
                            best_score = a[1]["AS"]

                    # Determine whether it is multimapping
                    alignment_paths = set()
                    for a in read_alignment:
                        if a[1]["AS"] == best_score:
                            # print(a[0][5])
                            alignment_paths.add(a[0][5])

                    if len(alignment_paths) > 1:
                        best_alignment[1]["at"] = "10"
                    else:
                        best_alignment[1]["at"] = "00"

                    same_read_pair[ri] = best_alignment

                # print(same_read_pair[0])
                # print(same_read_pair[1])
                # print(same_read_pair[2])
                # print()

                # print(x)
                yield (same_read_pair[1], same_read_pair[2])

                same_read_pair = [query_name, [], []]
                if r1flag:
                    same_read_pair[1].append(aln)
                if r2flag:
                    same_read_pair[2].append(aln)

        # TODO yield last read pair (edge case. Ugh)
        # return

    def get_alignments_parsed(self):

        for aln in self.get_alignments():
            for read_aln in aln:

                if read_aln is None:
                    continue

                regular_fields, tags = read_aln

                for i in [1, 2, 3, 6, 7, 8, 9, 10, 11]:
                    regular_fields[i] = int(regular_fields[i])

                regular_fields[5] = self.parse_path(regular_fields[5])

                if "cs" in tags:
                    tags["cs"] = self.parse_cs_tag(tags["cs"])

                # print(regular_fields)
                # print(tags)
                # print()
            yield aln


def phasing_main(gfa_fp, gaf_fp, maternal_walk, paternal_walk, output_fp, ignore_sex_chr=False):
    parental_strs = ["maternal", "paternal"]

    gfa_instance = GraphicalFragmentAssemblySegmentForPhasing()
    gfa_instance.parse(gfa_fp, maternal_walk, paternal_walk, ignore_sex_chr=ignore_sex_chr)

    gaf_instance = GraphAlignmentFormat(gaf_fp)

    output_fh = open(output_fp, "w")

    count_tmp = 0
    count_tmp1, count_tmp2 = 0, 0

    countm, countp, count_aln, countall = 0, 0, 0, 0
    for aln in gaf_instance.get_alignments_parsed():

        countall += 1

        # (read_index, M or P, supporting segment, position in segment)
        score_detail = []
        for read_index, read_aln in enumerate(aln):
            # aln is the read pair alignment
            # read_aln is the alignment for a single read

            if read_aln is None:
                continue
            count_aln += 1

            regular_fields, tags = read_aln
            path = regular_fields[5]
            # print(path)

            mat_unique_seg = set()
            pat_unique_seg = set()
            for seg, ori in path:
                if seg in gfa_instance._maternal_unique_segments:
                    mat_unique_seg.add(seg)
                if seg in gfa_instance._paternal_unique_segments:
                    pat_unique_seg.add(seg)

            # No unique segments, skip
            if len(mat_unique_seg) + len(pat_unique_seg) == 0:
                continue

            # Evaluate link orientation (if possible), if not, skip
            for parental_str in parental_strs:
                unique_seg = mat_unique_seg if parental_str == "maternal" else pat_unique_seg
                gfa_parental_info = gfa_instance._maternal_unique_segments if parental_str == "maternal" else gfa_instance._paternal_unique_segments

                for seg in list(unique_seg):
                    parent_seg_info, child_seg_info = gfa_parental_info[seg]

                    # Well, you can't check if there isn't any info
                    if parent_seg_info is None and child_seg_info is None:
                        continue

                    unique_seg_pos = None
                    for i, (this_seg_id, this_ori) in enumerate(path):
                        if this_seg_id == seg:
                            unique_seg_pos = i
                            break

                    forward_check = True
                    backward_check = True

                    upstream_segment_info_in_path = None
                    this_segment_info_in_path = path[unique_seg_pos]
                    downstream_segment_info_in_path = None

                    if unique_seg_pos > 0:
                        upstream_segment_info_in_path = path[unique_seg_pos - 1]
                    if unique_seg_pos < len(path) - 1:
                        downstream_segment_info_in_path = path[unique_seg_pos + 1]

                    # Check forward
                    if upstream_segment_info_in_path is not None and parent_seg_info is not None:
                        alignment_orientation = upstream_segment_info_in_path[1] * this_segment_info_in_path[1]

                        if upstream_segment_info_in_path[0] != parent_seg_info[0]:
                            forward_check = False

                            if alignment_orientation != parent_seg_info[1]:
                                forward_check = False

                    if downstream_segment_info_in_path is not None and child_seg_info is not None:
                        alignment_orientation = this_segment_info_in_path[1] * downstream_segment_info_in_path[1]

                        if downstream_segment_info_in_path[0] != child_seg_info[0]:
                            backward_check = False

                            if alignment_orientation != child_seg_info[1]:
                                backward_check = False

                    # Check backward
                    if downstream_segment_info_in_path is not None and parent_seg_info is not None:
                        alignment_orientation = this_segment_info_in_path[1] * downstream_segment_info_in_path[1]

                        if downstream_segment_info_in_path[0] != parent_seg_info[0]:
                            backward_check = False

                            if alignment_orientation != parent_seg_info[1]:
                                backward_check = False

                    if upstream_segment_info_in_path is not None and child_seg_info is not None:
                        alignment_orientation = upstream_segment_info_in_path[1] * this_segment_info_in_path[1]

                        if upstream_segment_info_in_path[0] != child_seg_info[0]:
                            backward_check = False

                            if alignment_orientation != child_seg_info[1]:
                                backward_check = False

                    if not (forward_check or backward_check):
                        # It means altough it aligns to a unique segment, the "context" isn't right
                        unique_seg.remove(seg)

            # No unique segment passed link check, skip
            if len(mat_unique_seg) + len(pat_unique_seg) == 0:
                continue

            # Evaluate alignment. If the read was aligned to a unique segment, just making sure the read is not mismatches
            # Rebuild the alignment block

            alignment_path_length = read_aln[0][6]
            alignment_path_length_reconstructed = 0
            for seg, ori in path:
                alignment_path_length_reconstructed += gfa_instance.get_sequence_length_by_segment_ID(seg)

            count_tmp1 += 1
            # print(alignment_path_length, alignment_path_length_reconstructed)
            if alignment_path_length != alignment_path_length_reconstructed:
                # This is known issue with VG giraffe. The alignment path length may not be correct. Skip this read when it happens.
                """
                print("Alignment path length mismatch")
                print(read_aln[0][5])
                print(alignment_path_length, alignment_path_length_reconstructed)
                for seg, ori in path:
                    print(seg, ori, gfa_instance.get_sequence_length_by_segment_ID(seg))
                print()
                print(alignment_path_length)
                print()
                count_tmp2 += 1"""
                count_tmp2 += 1
                # print(count_tmp2 / count_tmp1)
                continue
                # raise RuntimeError

            pos_in_alignment_block = 0
            pos_in_alignment_path = read_aln[0][7]
            pos_in_read = read_aln[0][2]
            # print(pos_in_alignment_block, pos_in_alignment_path, pos_in_read)

            current_segment_index = 0
            current_segment = path[current_segment_index][0]
            current_segment_orientation = path[current_segment_index][1]
            current_segment_length = gfa_instance.get_sequence_length_by_segment_ID(current_segment)
            pos_in_segment = read_aln[0][7]
            if current_segment_orientation == -1:
                pos_in_segment = current_segment_length - pos_in_segment - 1

            try:
                # Why does it happen?
                # When the read soft-clipped, and the true alignment of the starting segment is too short, it might extend to the previous one.
                assert pos_in_segment >= 0
                assert pos_in_segment < current_segment_length
            except:
                continue
                print(regular_fields)
                for seg, ori in path:
                    print(seg, ori, gfa_instance.get_sequence_length_by_segment_ID(seg))

            # JUST FOR DEBUGGING
            """ref_seq = ""
            read_seq = ""
            alignment_seq = ""
            current_segment_seq = gfa_ss_instance.get_sequence_by_segment_ID(current_segment)

            read_seq_original = read_data._data[read_aln[0][0]][read_index]
            # STOP DEBUGGING """

            alignment_block = read_aln[1]["cs"]
            # print(alignment_block)
            # print(read_aln[0][5])
            alignment_block_longer_than_path = False
            for atype, atype_length in alignment_block:
                # How do we deal with each alignment type
                # M: match. All good, add some score
                # X: mismatch. If it is aligned to a unique segment, then not add score for the mismatch.
                # I: insertion. If it is aligned to a unique segment, then not add score for the insertion portion.
                # D: deletion. If it is aligned to a unique segment, then not add score for the deletion portion for sure.

                if atype == "X":
                    atype_length = 1
                for ai in range(atype_length):

                    if alignment_block_longer_than_path:
                        break

                    if atype in "MX":

                        if atype == "M":
                            is_mat_unique_seg = current_segment in mat_unique_seg
                            is_pat_unique_seg = current_segment in pat_unique_seg
                            if is_mat_unique_seg:
                                score_detail.append((read_index, "M", current_segment, pos_in_segment))
                            if is_pat_unique_seg:
                                score_detail.append((read_index, "P", current_segment, pos_in_segment))

                        """
                        if atype == "M":
                            alignment_seq += " "
                        else:
                            alignment_seq += "X"
                        
                        read_seq += read_seq_original[pos_in_read]

                        # Hmmm
                        # print(pos_in_alignment_block)
                        ref_nucleotide = current_segment_seq[pos_in_segment]
                        if current_segment_orientation == -1:
                            ref_nucleotide = nucleotide_complement(ref_nucleotide)
                        ref_seq += ref_nucleotide
                        """

                        # Hmmm orientation of the path handling. Annoying
                        pos_in_segment += current_segment_orientation

                        next_segment_flag = False
                        if current_segment_orientation == 1 and pos_in_segment == current_segment_length:
                            next_segment_flag = True
                        elif current_segment_orientation == -1 and pos_in_segment == -1:
                            next_segment_flag = True
                        else:
                            # print(pos_in_segment, current_segment_length, current_segment_orientation)
                            pass

                        if next_segment_flag:

                            current_segment_index += 1
                            if current_segment_index == len(path):
                                # It could happen when the right end of the read aligns to the end of the path.
                                # Or VG report is not correct.
                                alignment_block_longer_than_path = True
                                continue
                                pass

                            current_segment = path[current_segment_index][0]
                            current_segment_orientation = path[current_segment_index][1]
                            current_segment_length = gfa_instance.get_sequence_length_by_segment_ID(current_segment)
                            # current_segment_seq = gfa_ss_instance.get_sequence_by_segment_ID(current_segment)
                            pos_in_segment = 0
                            if current_segment_orientation == -1:
                                pos_in_segment = current_segment_length - 1

                        pos_in_alignment_path += 1
                        pos_in_read += 1
                        pos_in_alignment_block += 1




                    elif atype == "I":

                        # alignment_seq += "+"
                        # read_seq += read_seq_original[pos_in_read]
                        # ref_seq += " "

                        # pos_in_alignment_path += 1
                        pos_in_read += 1
                        pos_in_alignment_block += 1





                    elif atype == "D":
                        """
                        alignment_seq += "-"
                        read_seq += " "

                        # Hmmm
                        ref_nucleotide = current_segment_seq[pos_in_segment]
                        if current_segment_orientation == -1:
                            ref_nucleotide = nucleotide_complement(ref_nucleotide)
                        ref_seq += ref_nucleotide
                        """

                        # Hmmm orientation of the path handling. Annoying
                        pos_in_segment += current_segment_orientation

                        next_segment_flag = False
                        if current_segment_orientation == 1 and pos_in_segment == current_segment_length:
                            next_segment_flag = True
                        elif current_segment_orientation == -1 and pos_in_segment == -1:
                            next_segment_flag = True
                        else:
                            # print(pos_in_segment, current_segment_length, current_segment_orientation)
                            pass

                        if next_segment_flag:

                            current_segment_index += 1
                            if current_segment_index == len(path):
                                alignment_block_longer_than_path = True
                                continue
                                pass

                            current_segment = path[current_segment_index][0]
                            current_segment_orientation = path[current_segment_index][1]
                            current_segment_length = gfa_instance.get_sequence_length_by_segment_ID(current_segment)
                            # current_segment_seq = gfa_ss_instance.get_sequence_by_segment_ID(current_segment)
                            pos_in_segment = 0
                            if current_segment_orientation == -1:
                                pos_in_segment = current_segment_length - 1

                        pos_in_alignment_path += 1
                        # pos_in_read += 1
                        pos_in_alignment_block += 1



        # Make the decision here
        # Exam the score_detail, if the read pair can be phased, output the phasing information
        # 1. Output the phasing information: NA, maternal, paternal, most likely maternal, most likely paternal, ambiguous
        # 2. Output unique mapped or multimapped
        # if len(score_detail) == 0:
        #    continue
        mat_score, pat_score = 0, 0
        score_sort = {}

        for read_index, parental, seg, pos in score_detail:
            if seg not in score_sort:
                score_sort[seg] = {}
            if pos not in score_sort[seg]:
                score_sort[seg][pos] = [None, 0, 0]
            if score_sort[seg][pos][0] is None:
                score_sort[seg][pos][0] = parental
            score_sort[seg][pos][read_index + 1] = 1

        for seg in score_sort.keys():
            for pos, (parental, read1score, read2score) in score_sort[seg].items():
                s = max(read1score, read2score)
                if parental == "M":
                    mat_score += s
                if parental == "P":
                    pat_score += s

        read_class = "NA"
        class_confidence = None
        if mat_score == 0 and pat_score == 0:
            read_class = "NA"
        elif mat_score > pat_score:
            read_class = "maternal"
            class_confidence = (mat_score - pat_score, pat_score)
            countm += 1
        elif pat_score > mat_score:
            read_class = "paternal"
            class_confidence = (pat_score - mat_score, mat_score)
            countp += 1
        else:
            read_class = "ambiguous"

        read_name = "*"
        cigar1 = "*"
        cigar2 = "*"
        mp1 = "*"
        mp2 = "*"
        for read_index, read_aln in enumerate(aln):
            if read_aln is None:
                continue
            read_name = read_aln[0][0]

            cigarx = read_aln[1]["cs"]
            mpx = read_aln[1]["at"]

            cigar_str = ""
            for atype, atype_length in cigarx:
                if atype == "X":
                    atype_length = 1

                cigar_str += f"{atype_length}{atype}"

            mp = "U"
            if mpx[0] == "1":
                mp = "M"

            if read_index == 0:
                cigar1 = cigar_str
                mp1 = mp
            else:
                cigar2 = cigar_str
                mp2 = mp

        line = [read_name, read_class, mat_score, pat_score, mp1, mp2, cigar1, cigar2]
        line_str = "\t".join([str(x) for x in line]) + "\n"
        output_fh.write(line_str)
        # print(line)
        # print(score_detail)
        # print()

        # print(score_detail)

        count_tmp += 1
        if count_tmp > 100:
            pass  # break

        if countall % 10000 == 0:
            pass
            # print(countm, countp, count_aln, countall)

    output_fh.close()
    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Genome graph based phasing tool')
    # Required arguments
    # path to genome graph in GFA format
    # path to alignment file in GAF format
    # Maternal walk
    # Paternal walk
    # output file name

    # Optional arguments
    # verbose
    # ignore sex chromosomes # true or false

    # Required arguments
    parser.add_argument('-gfa', type=str, help='Path to genome graph in GFA format', required=True)
    parser.add_argument('-gaf', type=str, help='Path to alignment file in GAF format', required=True)
    parser.add_argument('-maternal', type=str, help='Maternal walk', required=True)
    parser.add_argument('-paternal', type=str, help='Paternal walk', required=True)
    parser.add_argument('-output', type=str, help='Output file name', required=True)


    # Optional arguments
    parser.add_argument('-verbose', help='Print verbose output')
    parser.add_argument('-ignore_sex_chr', help='Ignore the sex chromosomes', default="False")


    # Parse the arguments
    args = parser.parse_args()

    gfa_fp = args.gfa
    gaf_fp = args.gaf
    maternal_walk_str = args.maternal
    paternal_walk_str = args.paternal
    output_fp = args.output

    verbose = args.verbose
    ignore_sex_chr = args.ignore_sex_chr


    maternal_walk = maternal_walk_str.split(":")
    paternal_walk = paternal_walk_str.split(":")

    assert len(maternal_walk) == 2
    assert len(paternal_walk) == 2

    maternal_walk[1] = int(maternal_walk[1])
    paternal_walk[1] = int(paternal_walk[1])

    # print(maternal_walk, paternal_walk)


    # Fact check the arguments
    if not os.path.exists(gfa_fp):
        print(f"Genome graph file {gfa_fp} does not exist", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(gaf_fp):
        print(f"Alignment file {gaf_fp} does not exist", file=sys.stderr)
        sys.exit(1)

    if ignore_sex_chr.lower() in ["true", "t", "yes", "y"]:
        ignore_sex_chr = True
    elif ignore_sex_chr.lower() in ["false", "f", "no", "n"]:
        ignore_sex_chr = False
    else:
        print(f"Invalid value for ignore_sex_chr {ignore_sex_chr}", file=sys.stderr)

    # print(gfa_fp, gaf_fp, maternal_walk_str, paternal_walk_str, output_fp, ignore_sex_chr)

    phasing_main(gfa_fp, gaf_fp, maternal_walk, paternal_walk, output_fp, ignore_sex_chr=ignore_sex_chr)

































