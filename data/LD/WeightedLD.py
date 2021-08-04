#!/usr/bin/env python3
# Oscar Charles & Joseph Roberts 2021
# This code calculates sequence weights using the Henikoff formula from a multiple sequence alignment
# Usage: python WeightedLD.py --input alignment.fasta

import logging
import argparse
from pathlib import Path
from Bio import AlignIO
import numpy as np
import sys
import re


logging.basicConfig(
    format='[%(levelname)s] %(asctime)s %(message)s',
    level=logging.ERROR,
    datefmt='%Y-%m-%d %H:%M:%S'
)


def read_fasta(filename: Path) -> np.ndarray:
    """
    Converts a fasta DNA alignment file to a hashed integer matrix
    """
    raw_alignment = AlignIO.read(filename, "fasta")

    alignment_chars = np.zeros(
        (len(raw_alignment), raw_alignment.get_alignment_length()), dtype='<U1')

    # forces lowercase IUPAC characters
    for (iSeq, seq) in enumerate(raw_alignment):
        alignment_chars[iSeq, :] = list(str(seq.seq).lower())

    alignment = np.full_like(alignment_chars, fill_value=5, dtype=np.uint8)
    alignment[alignment_chars == 'a'] = 0
    alignment[alignment_chars == 'c'] = 1
    alignment[alignment_chars == 'g'] = 2
    alignment[alignment_chars == 't'] = 3
    alignment[alignment_chars == '-'] = 4

    return alignment


def compute_variable_sites(alignment: np.ndarray, min_acgt: float, min_variability: float) -> np.ndarray:
    """
    returns which sites have both sufficient data and a significant enough quantity of a minor symbol

    Args:
        alignment: 2D numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (A, C, G, T, -, ambiguous) is
              represented by (0, 1, 2, 3, 4, 5)
        min_acgt: The minimum fraction of sequences which must have A/C/G/T
            symbols at a given site
        min_variability: The minimum fraction of sequences which have the minor
            symbol at a given site, ignoring sequences that have missing or
            ambiguous symbols.

    Returns:
        2 member tuple of 1D boolean arrays of length n_sites, which is True for each site that
        meets the criteria.
        1 relevant for Henikoff calculations
        2 relevant for LD calculations
    """

    # For each site, the fraction of sequences which have a concrete symbol at that position
    concrete_fraction = (alignment < 4).sum(axis=0) / alignment.shape[0]

    # Whether a given site has enough data to be considered for further processing
    sufficient_data = concrete_fraction > min_acgt

    # For each site, how many sequences have a given symbol at that site
    acgt_counts = np.array([(alignment == x).sum(axis=0)
                           for x in (0, 1, 2, 3, 4)])

    # Whether each site has multiple non-ambiguous symbols - deprecated
    # multiple_non_ambiguous = np.sum(acgt_counts > 0, axis=0) > 1
    major_counts = acgt_counts.max(axis=0)
    minor_counts = acgt_counts.sum(axis=0) - major_counts

    # sites with any variation
    multiple_non_ambiguous = minor_counts > 0
    minor_fraction = np.zeros(alignment.shape[1])
    minor_fraction[multiple_non_ambiguous] = minor_counts[multiple_non_ambiguous] / \
        (major_counts[multiple_non_ambiguous] +
         minor_counts[multiple_non_ambiguous])

    # Does the minor symbol occur at a high enough frequency
    has_min_variability = minor_fraction >= min_variability

    # invariant sites should be included, so we should only exclude sites with poor coverage
    return_hk_varsites = sufficient_data

    # has enough variability to return useful LD data
    return_ld_varsites = sufficient_data & has_min_variability

    return return_hk_varsites, return_ld_varsites


def henikoff_weighting(alignment: np.ndarray) -> np.ndarray:
    """
    The Henikoff weighting for each sequence in the given alignment array.
    ref:
        Steven Henikoff and Jorja G. Henikoff (1994) "Position-based Sequence Weights"

    Args:
        alignment: 2D numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (A, C, G, T, -, ambiguous) is
              represented by (0, 1, 2, 3, 4, 5)
            - Only sites of interest are included

    Returns:
        A 1D numpy array of length n_seqs, containing the Henikoff weight of each sequence.
    """

    n_sites = alignment.shape[1]
    n_seqs = alignment.shape[0]

    # Mask for non-ambiguous bases (False where the base is ambiguous)
    ok_base = (alignment != 5)

    # For each site, the count of sequences that have a given symbol at that site
    # Eg count_base[1, 1234] contains the count of sequences that have symbol 1 at site 1234
    count_base = np.zeros((6, n_sites))
    for base in range(6):
        count_base[base, :] = (alignment == base).sum(axis=0)

    # For each site, the count of unique bases at that site
    unique_base = count_base[:5, :] > 0
    unique_base = np.sum(unique_base, axis=0)

    # The Henikoff weighting per site, per ok base
    site_contribution = np.zeros(alignment.shape)
    site_contribution[ok_base] = 1 / \
        (unique_base * count_base[alignment, np.arange(n_sites)])[ok_base]

    # For each ambiguous base, fill in the average contribution across all
    # sequences with concrete bases for that site
    site_contribution[~ok_base] = 0
    site_average_weight = site_contribution.sum(
        axis=0) / np.sum(count_base[:5, :], axis=0)
    site_contribution[~ok_base] = np.broadcast_to(
        site_average_weight, site_contribution.shape)[~ok_base]

    # For each sequence, the sum of contributions from each site
    weights = site_contribution.sum(axis=1)

    # Normalize such that the largest weight is exactly 1
    return weights / weights.max()


def ld(alignment, weights, site_map, r2_threshold=0.1):
    """
    Generate LD metrics; D, D', and R2

    Args:
        alignment: 2D numpy array where:
            - The first axis represents sequences
            - The second axis represents sites
            - Each element is an integer where (A, C, G, T, -, ambiguous) is
              represented by (0, 1, 2, 3, 4, 5)
            - Only sites of interest should be included
        weights: 1D numpy array of length n_seqs, containing the relative weight of each sequence
        site_map: Some object that can map site indices in the alignment argument to meaningful site indices

    Returns:
        Nothing, just prints the results to stdout in a tab-separated format
    """

    n_seqs = alignment.shape[0]
    n_sites = alignment.shape[1]

    # stdout headers
    print("site_a\tsite_b\tD\tD'\tr2")
    for first_site in range(n_sites - 1):
        logging.info("    Outer loop: %s/%s", first_site, n_sites)
        for second_site in range(first_site + 1, n_sites):
            # Form an array which contains all the sequences, but only the two target sites
            target_sites = alignment[:, (first_site, second_site)]
            # Remove all sequences with a bad symbol at either target site
            good_sequences = (target_sites < 5).all(axis=1)
            target_sites = target_sites[good_sequences, :]
            target_weights = weights[good_sequences]
            target_seqs = target_sites.shape[0]

            # For each of the sequence, the "major allele" is the symbol that occurs most frequently, dominantMinor the 2nd most
            # Whether the given sequence is equal to the major symbol at the given site
            target_sites_major = np.zeros_like(target_sites, dtype=np.bool8)
            target_sites_domMinor = np.zeros_like(target_sites, dtype=np.bool8)

            skip_site = False
            for site in (0, 1):
                unique_elements, counts = np.unique(
                    target_sites[:, site], return_counts=True)
                if len(unique_elements) <= 1:
                    # After removing the bad sequences, one or both of the
                    # sites may no longer be variable. Stop calculations here
                    # if that is the case.
                    skip_site = True
                # identifies positions with major_allele
                major_symbol = unique_elements[counts.argmax()]
                major_symbol = unique_elements[np.argsort(-counts)[0]]
                target_sites_major[target_sites[:, site]
                                   == major_symbol, site] = True
                if not skip_site:
                    # identifies positions with domMino_allele, if two are equal takes first
                    domMinor_symbol = unique_elements[np.argsort(-counts)[1]]
                    target_sites_domMinor[target_sites[:, site]
                                          == domMinor_symbol, site] = True
            if skip_site:
                continue

            # second round of filtering - remove anything thats not Major or domMinor
            # identify seqs with Maj or domMinor in both sites -> keep
            keep_first = target_sites_major[:, 0] + target_sites_domMinor[:, 0]
            keep_second = target_sites_major[:,
                                             1] + target_sites_domMinor[:, 1]
            keep = np.array([keep_first, keep_second]).all(axis=0)
            # filter again
            target_sites = target_sites[keep, :]
            # print(target_sites)
            target_weights = target_weights[keep]
            target_sites_major = target_sites_major[keep, ]
            target_seqs = target_sites.shape[0]
            total_weight = target_weights.sum()

            # ----- Calculate allele frequencies
            PA, PB = np.ma.masked_array(target_weights.reshape(-1, 1).repeat(
                2, axis=1), ~target_sites_major).sum(axis=0) / total_weight
            Pa, Pb = np.ma.masked_array(target_weights.reshape(-1, 1).repeat(
                2, axis=1), target_sites_major).sum(axis=0) / total_weight

            # after removing sequences which not Maj or dMin in both sites, we may want to skip site if a site is invariant
            if round(PA, 1) == 1.0:
                continue
            if round(PB, 1) == 1.0:
                continue

            # ----- predicted haplotype frequencies
            # When haplotype frequencies are equal to the product of their corresponding allele frequencies, then loci are in linkage equilibrium
            PAB = PA * PB
            PAb = PA * Pb
            PaB = Pa * PB
            Pab = Pa * Pb

            # ----- observed haplotype frequencies
            ld_obs = np.zeros(4)
            ld_obs[0] = target_weights[~target_sites_major[:, 0]
                                       & ~target_sites_major[:, 1]].sum()
            ld_obs[3] = target_weights[target_sites_major[:, 0]
                                       & target_sites_major[:, 1]].sum()
            ld_obs[1] = target_weights[~target_sites_major[:, 0]
                                       & target_sites_major[:, 1]].sum()
            ld_obs[2] = target_weights[target_sites_major[:, 0]
                                       & ~target_sites_major[:, 1]].sum()
            ld_obs = ld_obs / total_weight

            # ----- Caclulate D [Linkage Disequilibrium]
            # the vector is now as in the hahn molpopgen book and can be used to generate D values
            # ld_obs is pAB, p
            tD = np.zeros(4)
            tD[0] = PAB - ld_obs[3]  # MajMaj
            tD[1] = Pab - ld_obs[0]  # minmin
            tD[2] = -1 * (PAb - ld_obs[2])  # Majmin
            tD[3] = -1 * (PaB - ld_obs[1])  # minMaj
            D = (tD[0] + tD[1] + tD[2] + tD[3]) / \
                4  # they should be the same anyhow

            # normalised D = D'
            if D < 0:
                denominator = max([-PAb, -PaB])
            else:
                denominator = min([PAB, Pab])
            DPrime = D / denominator
            # calculate R2
            R2 = D**2 / (PA * Pa * PB * Pb)

            # --r2-threshold
            if(R2 < r2_threshold):
                continue

            # cat output
            print(
                f"{site_map[first_site]}\t{site_map[second_site]}\t{round(D, 4)}\t{round(DPrime, 4)}\t{round(R2, 4)}")


def handle_fasta(args):
    logging.info("Reading FASTA data from %s", args.input)
    # convert fasta character alignment to integer matrix
    alignment = read_fasta(args.input)
    logging.info(
        "Finished reading data. Sequence count: %s, Sequence length %s", *alignment.shape)

    logging.info("Computing sites of iterest (min_acgt=%s, min_variability=%s)",
                 args.min_acgt, args.min_variability)

    # identify informative sites
    var_sites_HK, var_sites_LD = compute_variable_sites(
        alignment, args.min_acgt, args.min_variability)
    logging.info("Found %s sites of interest", var_sites_LD.sum())

    # calculate Henikoff weights, do this each time its fast and logicflow otherwise is a bit complex
    weights = henikoff_weighting(alignment[:, var_sites_HK])

    # Trim down the alignment array to only include the sites of interest
    alignment = alignment[:, var_sites_LD]

    # Maps site indices in the trimmed down array to site indices in the original alignment
    site_map = np.where(var_sites_LD)[0]

    return alignment, site_map, weights


def handle_vcf(filename):
    # extract lines -> list
    logging.info("Reading VCF data from %s", filename)
    with open(filename, 'r') as f:
        data = f.read().split("\n")
    endline = len(data)

    # is there a header block? if yes return header line and data block
    is_vcf = False
    for i, line in enumerate(data):
        if re.search("#CHROM", line):
            is_vcf = True
            header = data[i]
            data = data[(i+1):endline]
            break
    if is_vcf == False:
        print(
            "No #CHROM header block identified")
        sys.exit(1)

    # is there enough dat to be meaningful?
    line = data[0].split("\t")
    if len(line) <= 12:
        print(
            "The VCF data contains too small a population, are you sure this is a multi VCF?")
        sys.exit(1)

    # is the data haploid or diploid?
    if(type(re.search(r"[0-2]|[0-2]", data[0])) == "NoneType"):
        # the vcf is haploid
        print(
            "Well this is awkward, we haven't implemented a haploid VCF reader yet")
        sys.exit(1)
    else:

        # the vcf is diploid
        # now split any diploid -> haploid
        for i, line in enumerate(data):
            # remove | delimiters in extra cols
            t = line.replace("|||", "")
            t = t.replace("||", "")
            t = re.sub(r"[^0-9]\|[^0-9]", "", t)
            t = re.sub(r"[^0-9]\|[^0-9]", "", t)
            # if any diplod calls are unphased (i.e. we do not know haplotype) treat as missing.
            t = re.sub(r"./.", ".|.", t)
            t = t.replace("|", "\t")
            # replace .  with proper missing flag
            t = re.sub(r"\.{1}", "4", t)

            t = t.split("\t")
            del t[2:9]
            del t[0]
            data[i] = t

        del data[len(data)-1]  # it is standard to leave the last line blank

    # we now have for each row the pos and haploid calls
    # extract site_map - large numbers
    site_map = np.array([item[0] for item in data], dtype=np.int64)

    # force everything to uint8 - remove site_map col
    alignment = np.array(data, dtype=np.uint8)
    alignment = np.delete(alignment, 0, axis=1)
    # so that its the same format as an alignment - for compatability with other functions
    alignment = np.rot90(alignment)

    # calculate Henikoff weights, do this each time its fast and logicflow otherwise is a bit complex
    weights = henikoff_weighting(alignment)

    logging.info(
        "Finished reading data. Sequence count: %s, Sequence length %s", *alignment.shape)
    return alignment, site_map, weights


def main(args):
    filename = str(args.input)
    if filename.endswith('.vcf'):
        # with a VCF all sites should have been filtered to be variant, less Henikoff information, but will have to do
        alignment, site_map, weights = handle_vcf(filename)
    else:
        alignment, site_map, weights = handle_fasta(args)

    # default behaviour is Henikoff weighting, can be disabled
    if args.unweighted:
        logging.info("Unweighted")
        weights = np.zeros(alignment.shape[0], dtype=np.uint8)
        weights[weights == 0] = 1
    else:
        logging.info("Computing Henikoff weights for each sequence")
        # print Henikoff weights to table
        if args.weights_output != None:
            np.savetxt(args.weights_output, weights, delimiter='\t')

    logging.info("Computing the LD parameters")

    ld(alignment, weights, site_map, args.r2_threshold)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WeightedLD computation tool")
    parser.add_argument("--input", type=Path,
                        help="A multiple sequence alignment in FASTA format, or multi sample VCF", required=True)
    parser.add_argument("--min-acgt", type=float, default=0.8,
                        help="Sets a minimum fraction of A,C,G & T required for a site to be considered in LD and \
                            weighting calculations. Increase to account for poor sequence coverage.")
    parser.add_argument("--min-variability", type=float, default=0.02,
                        help="The minimum (dominant) minor allele fraction for a site to be considered in LD calculations")
    parser.add_argument("--r2-threshold", type=float, default=0.1,
                        help="Minimum value of R2 for a pairwise site comparion to be included in the output")
    parser.add_argument("--weights-output", type=Path, required=False,
                        help="Filename to write the per-sequence weights to, in Tab Separated Value format")
    parser.add_argument("--unweighted", action='store_true', default=False,
                        help="Use unit weights instead of Henikoff weights")

    args = parser.parse_args()
    main(args)
