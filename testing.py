import euler

"""Assembles a set of kmers using the Eulerian path algorithm.
    Args:
        kmers: a set of kmers (DNA strings with length = k)
    Returns:
        A shortest superstring with k-mer spectrum equal to the input set of kmers
"""
def euler_assemble(kmers):
    return euler.euler_assembly(kmers)

# TEST: euler_assemble returns a string
sanity_test_kmers = {'AC', 'AG', 'AT', 'CA', 'CC',
                     'CT', 'GA', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'}
assert isinstance(euler_assemble(sanity_test_kmers), str)
print("SUCCESS: euler_assemble returns a string passed!")

# TEST: euler_assemble returns a superstring of the kmers
def is_superstring(s, kmers):
    return all(kmer in s for kmer in kmers)
assert is_superstring(euler_assemble(sanity_test_kmers), sanity_test_kmers)
print("SUCCESS: euler_assemble returns a superstring passed!")

# TEST: euler_assemble returns a shortest string with spectrum identical to the set of input kmers
def spectrum(s, k):
    return {s[i: i + k] for i in range(len(s) - k + 1)}
def has_spectrum(s, k, kmers):
    assert isinstance(s, str)
    return spectrum(s, k) == kmers and len(kmers) == (len(s) - k + 1)
assert has_spectrum(euler_assemble(sanity_test_kmers), 2, sanity_test_kmers)
print("SUCCESS: euler_assemble returns a string with a matching spectrum!")

# TEST: euler_assemble_3_3
kmers_3_3 = {"GGA"}
assert has_spectrum(euler_assemble(kmers_3_3), 3, kmers_3_3)
print("SUCCESS: euler_assemble_3_3 passed!")

# TEST: euler_assemble_3_4
kmers_3_4 = {"AGT", "GTG"}
assert has_spectrum(euler_assemble(kmers_3_4), 3, kmers_3_4)
print("SUCCESS: euler_assemble_3_4 passed!")

# TEST: euler_assemble_3_5
kmers_3_5 = {"CGC", "GCT", "TCG"}
assert has_spectrum(euler_assemble(kmers_3_5), 3, kmers_3_5)
print("SUCCESS: euler_assemble_3_5 passed!")

# TEST: euler_assemble_3_10
kmers_3_10 = {"ATG", "CGC", "CGT", "GAT", "GCG", "GTG", "TGA", "TGC"}
assert has_spectrum(euler_assemble(kmers_3_10), 3, kmers_3_10)
print("SUCCESS: euler_assemble_3_10 passed!")

# TEST: euler_assemble_3_20
kmers_3_20 = {'CTT', 'AGC', 'GCA', 'TGC', 'GGC', 'CAG', 'CCA', 'GAT', 'ATG',
              'TCT', 'GCG', 'CAT', 'TCC', 'TTT', 'CGA', 'GCT', 'ATC', 'CTC'}
assert has_spectrum(euler_assemble(kmers_3_20), 3, kmers_3_20)
print("SUCCESS: euler_assemble_3_20 passed!")

# TEST: euler_assemble_hidden_1
kmers_3_60 = {'CCC', 'AGT', 'GGG', 'GTC', 'GTT', 'TCG', 'TGC', 'GAA', 'CGG', 'TCT',
              'GCG', 'CAT', 'CAC', 'ATC', 'CTT', 'TTC', 'AGG', 'CAG', 'CCA', 'CCG',
              'TAC', 'GAC', 'TGA', 'TTA', 'TTT', 'ACG', 'TTG', 'CGT', 'ACC', 'CTC',
              'CTG', 'GTA', 'GCA', 'TAT', 'AAG', 'GGT', 'ATG', 'TCC', 'ATA', 'CTA',
              'ATT', 'GCT', 'AAT', 'AGC', 'ACT', 'CCT', 'CGC', 'GGC', 'AGA', 'TCA',
              'TGG', 'TAG', 'GCC', 'GGA', 'TAA', 'GAG', 'CGA', 'GAT'}
assert has_spectrum(euler_assemble(kmers_3_60), 3, kmers_3_60)
print("SUCCESS: euler_assemble_hidden_1 passed!")

# TEST: euler_assemble_hidden_2
kmers_4_60 = {'ACGC', 'CAAG', 'GGCT', 'GTGG', 'CACG', 'GGCA', 'CTTA', 'AACG', 'GTTT', 'TCAA',
              'CTCT', 'ACAC', 'TTTT', 'TAAA', 'TGGC', 'GCAT', 'GGTG', 'AAAC', 'TAAC', 'CAGG',
              'TTAA', 'ACTG', 'ACAG', 'CGCG', 'GGAC', 'TTCA', 'TGCT', 'TGGG', 'TTTC', 'TCGT',
              'AACA', 'ATCA', 'GTGC', 'AAGG', 'TTCG', 'ATTC', 'TATT', 'GACT', 'GCGT', 'GTAT',
              'CGTA', 'TCTA', 'AGGC', 'CGGT', 'GCTC', 'CTAA', 'ACGG', 'TCAC', 'AGGA', 'GCTT',
              'CATC', 'CGTG', 'AGCA', 'CTGG', 'CACA', 'GGTT', 'GGGT'}
assert has_spectrum(euler_assemble(kmers_4_60), 4, kmers_4_60)
print("SUCCESS: euler_assemble_hidden_2 passed!")

# TEST: euler_assemble_hidden_3
kmers_5_100 = {'TGATC', 'GAACA', 'GACAC', 'GACTG', 'CCCGC', 'GAGAC', 'GACTT', 'CCAAC', 'GTATC', 'GGAAC',
               'CGACT', 'TAGGA', 'GCGAC', 'TGCAG', 'TTCCC', 'TTAGG', 'ACGAC', 'AACTG', 'TGAGT', 'CTGCC',
               'ATCGG', 'GAGTG', 'GAACT', 'GTGAT', 'TAGTG', 'ATACG', 'TCGGT', 'CATAG', 'ATCTG', 'GGTTA',
               'TCGCG', 'CCGCG', 'ACTGC', 'ATAGT', 'CGGAA', 'TCTGA', 'AACGC', 'CGGGA', 'ACTTC', 'CTGAG',
               'ATCGA', 'AATCG', 'GCGTA', 'GGACT', 'ACCGG', 'ATACT', 'CCCCG', 'AATAC', 'CCGGA', 'CGGTT',
               'CATAC', 'ACACC', 'TGTGC', 'ACGCG', 'CGCGG', 'AACAT', 'ACTGA', 'ATCGC', 'AGGAC', 'CTTCC',
               'GAATA', 'AGCAT', 'GTGTG', 'TCCCC', 'GTGCA', 'CACCG', 'CGCGT', 'GGAGA', 'GGGAG', 'GCATA',
               'ACATC', 'TCGAA', 'GCAGC', 'GCGGG', 'GACAT', 'AGTGA', 'AGTGT', 'ACATA', 'CAACG', 'CTGAA',
               'CGTAT', 'CATCT', 'CGCGA', 'GTTAG', 'GCCAA', 'GGAAT', 'GATCG', 'CGACA', 'GAATC', 'AGACA',
               'CAGCA', 'TGCCA', 'TACGA', 'CGAAC', 'TATCG', 'TGAAT'}
assert has_spectrum(euler_assemble(kmers_5_100), 5, kmers_5_100)
print("SUCCESS: euler_assemble_hidden_3 passed!")

# TEST: euler_assemble_hidden_4
kmers_10_1000 = {line.strip() for line in open("10_1000.kmers.txt")}
assert has_spectrum(euler_assemble(kmers_10_1000), 10, kmers_10_1000)
print("SUCCESS: euler_assemble_hidden_4 passed!")

# TEST: euler_assemble_hidden_5
kmers_50_500 = {line.strip() for line in open("50_500.kmers.txt")}
assert has_spectrum(euler_assemble(kmers_50_500), 50, kmers_50_500)
print("SUCCESS: euler_assemble_hidden_5 passed!")
