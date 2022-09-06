def get_seq_comp(sequence):
    dna_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    comp_strand = ""
    for c in sequence:
        comp_strand += dna_dict[c]
    return comp_strand


def cut_original_plasmid_strand(sequence, rest_site):
    cut = sequence.find(rest_site)
    return sequence[:cut + 1], sequence[cut + 1:]


def cut_comp_plasmid_strand(comp_seq, comp_rest_site):
    cut = comp_seq.rfind(comp_rest_site)
    return comp_seq[0:cut + 5], comp_seq[cut + 5:]


def cut_GFP_strand(sequence, rest_site1, rest_site2):
    cut1 = sequence.find(rest_site1)
    cut2 = sequence.rfind(rest_site2)
    return sequence[cut1 + 1:cut2 + 1]


def get_GFP_sticky_ends(GFP_seq, L_rest_site, R_rest_site):
    GFP_seq_comp = get_seq_comp(GFP_seq)
    L_rest_site_comp = get_seq_comp(L_rest_site)
    R_rest_site_comp = get_seq_comp(R_rest_site)

    chain1 = cut_GFP_strand(GFP_seq, L_rest_site, R_rest_site)
    chain2 = cut_GFP_strand(GFP_seq_comp[::-1], R_rest_site_comp[::-1],
                            L_rest_site_comp[::-1])[::-1]
    return chain1, chain2


def plasmid_GFP_ligation(plasmid_strand, plasmid_rest_site, GFP_strand,
                         GFP_left_rest_Site, GFP_right_rest_site):
    plasmid_strand_part1, plasmid_strand_part2 = cut_original_plasmid_strand(
        plasmid_strand, plasmid_rest_site)
    comp_plasmid_strand = get_seq_comp(plasmid_strand)
    comp_plasmid_rest_site = get_seq_comp(plasmid_rest_site)
    comp_plasmid_part1, comp_plasmid_part2 = cut_comp_plasmid_strand(
        comp_plasmid_strand, comp_plasmid_rest_site)
    L_GFP_sticky_end, R_GFP_sticky_end = get_GFP_sticky_ends(GFP_strand,
                                                             GFP_left_rest_Site,
                                                             GFP_right_rest_site)
    strand1 = plasmid_strand_part1 + L_GFP_sticky_end + plasmid_strand_part2
    strand2 = comp_plasmid_part1 + R_GFP_sticky_end + comp_plasmid_part2

    return strand1, strand2


def process_input(file):
    data = list()
    with open(file, "r", encoding="utf-8") as input_file:
        for line in input_file:
            data.append(line.strip())
    return data


def run():
    file_name = input()
    sequences = process_input(file_name)
    plasmid_strand = sequences[0]
    plasmid_rest_site = sequences[1]
    GFP_strand = sequences[2]
    GFP_left_rest_site, GFP_right_rest_site = sequences[3].split()

    plasmid_strand1, plasmid_strand2 = plasmid_GFP_ligation(plasmid_strand,
                                                            plasmid_rest_site,
                                                            GFP_strand,
                                                            GFP_left_rest_site,
                                                            GFP_right_rest_site
                                                            )

    print(plasmid_strand1)
    print(plasmid_strand2)


if __name__ == '__main__':
    run()
