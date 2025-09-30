<file_path>
PhylogenicTree/docs/rds_inspection_summary.txt
</file_path>

<edit_description>
Update the RDS inspection summary with the actual output from running the script
</edit_description>

# RDS Inspection Summary

This summary contains the actual output from running the inspection script `scripts/phase_01_explore_validate/inspect_rds.R` on the RDS files in `data/interim/rds/`. The script loads each RDS file, prints a `str()` summary, and checks key fields with `head()`.

Rscript inspect_rds.R
=== Loading: pairing_result.rds ===
str() summary:
List of 52
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 $ :List of 1
 - attr(*, "dim")= int [1:2] 26 2
NULL

Key head() checks:

=== Loading: rec_cluster.rds ===
str() summary:
List of 2
 $ distm   : num [1:42190, 1:42190] 0 0.989 0.971 0.971 0.974 ...
 $ dataname: chr [1:42190] "GCF_000165835.1_AXYL_RS02075" "CRBC_G2531_ctg3_554" "CRBC_G2530_ctg6_182" "CRBC_G2627_ctg4_296" ...
NULL

Key head() checks:
head(dataname):  GCF_000165835.1_AXYL_RS02075, CRBC_G2531_ctg3_554, CRBC_G2530_ctg6_182, CRBC_G2627_ctg4_296, GCF_021496925.1_YS110_RS20765, CRBC_G5874_ctg310_26
dim(distm):  42190 x 42190

=== Loading: rec.rds ===
str() summary:
List of 10
 $ count       : num [1(1d)] 42190
 $ foldername  : chr [1:42190] "Burkholderia_cepacia_ATCC_39356" "Burkholderia_cepacia_ATCC_39356" "Burkholderia_cepacia_ATCC_39356" "Burkholderia_cepacia_ATCC_39356" ...
 $ fragmentname: chr [1:42190] "DEFINITION Burkholderia_cepacia_ATCC_39356." "DEFINITION Burkholderia_cepacia_ATCC_39356." "DEFINITION Burkholderia_cepacia_ATCC_39356." "DEFINITION Burkholderia_cepacia_ATCC_39356." ...
 $ boarders    : num [1:42190, 1:2] 1792443 3479322 4405967 4628948 4898646 ...
 $ seqs        : chr [1:42190] "MEWATGTRVRAIAAAASVAFGMAAGHAYAQTAPAVNAGAAASAGNVQNGATSGTLPAINVNAGSEGDGTVGLVAKRSRTGTKTDTSINEIPQTINVVTAQQIEMTGATDVN"| __truncated__ "MKETDGAKRRFAARRAGVTLCFGGMLGGWLTPASAQVERAGNVDAPTLAPIVVVGTTPLLGIGTPLSRVPANVQTIRGEDIMRQHGSVLTDYFEKNVSSVDINEAQGNPYQ"| __truncated__ "MLMETGAARRAWAVAAGGTLCVAAAGGAHAQAPNAGATTAVLPPIDVTGNAGGSSVGLVGLRTAAGTKTDTPVAEIPQTTNIVTAQQIEMTGAADLNQALRYVPGFATFGA"| __truncated__ "MPSRRPPRRARTLRPWRASLPALITLCVASGAHADAEPAATASPAAPSGSSTPAPPERELPTISVSASAATDPTVGYQPRTSSVAGGDDRPLKEIPQSVAVVSSSVMQDQQ"| __truncated__ ...
 $ tag         : chr [1:42190] "ctg1_1598" "ctg1_3141" "ctg1_3952" "ctg1_4119" ...
 $ recname     : chr [1:42190] "Burkholderia_cepacia_ATCC_39356_ctg1_1598" "Burkholderia_cepacia_ATCC_39356_ctg1_3141" "Burkholderia_cepacia_ATCC_39356_ctg1_3952" "Burkholderia_cepacia_ATCC_39356_ctg1_4119" ...
 $ domloc      : num [1:42190, 1:2] 88 67 74 92 82 70 61 42 79 54 ...
 $ domseq      : chr [1:42190] "GAAASAGNVQNGATSGTLPAINVNAGSEGDGTVGLVAKRSRTGTKTDTSINEIPQTINVVTAQQIEMTGATDVNAALRYVPGFSSYGSDNRSDWYAALRGFTPTAYVNGLQ"| __truncated__ "GVTLCFGGMLGGWLTPASAQVERAGNVDAPTLAPIVVVGTTPLLGIGTPLSRVPANVQTIRGEDIMRQHGSVLTDYFEKNVSSVDINEAQGNPYQTDVNYRGFTASPIVGT"| __truncated__ "AAGGAHAQAPNAGATTAVLPPIDVTGNAGGSSVGLVGLRTAAGTKTDTPVAEIPQTTNIVTAQQIEMTGAADLNQALRYVPGFATFGADSRTDWYAALRGFTPTLYVDGVP"| __truncated__ "ASPAAPSGSSTPAPPERELPTISVSASAATDPTVGYQPRTSSVAGGDDRPLKEIPQSVAVVSSSVMQDQQARSLDDVLGNISGVTQTNTLGGTRDAFIKRGFGSNNDGSVL"| __truncated__ ...
 $ group       : num [1:42190, 1] 4667 4112 4661 4371 1267 ...
NULL

Key head() checks:
head(recname):  Burkholderia_cepacia_ATCC_39356_ctg1_1598, Burkholderia_cepacia_ATCC_39356_ctg1_3141, Burkholderia_cepacia_ATCC_39356_ctg1_3952, Burkholderia_cepacia_ATCC_39356_ctg1_4119, Burkholderia_cepacia_ATCC_39356_ctg1_4363, Burkholderia_cepacia_ATCC_39356_ctg1_4511
head(group):  4667, 4112, 4661, 4371, 1267, 4956

=== Loading: syn_cluster.rds ===
str() summary:
List of 2
 $ distm   : num [1:2687, 1:2687] 0 0 0.001563 0.000781 0.000781 ...
 $ dataname: chr [1:2687] "GCF_000317955.1_NC_018829.1.region009" "GCF_002859225.1_NZ_CP018761.1.region009" "GCF_003186595.1_NZ_CP020645.1.region009" "GCF_008693665.1_NZ_CP044089.1.region005" ...
NULL

Key head() checks:
head(dataname):  GCF_000317955.1_NC_018829.1.region009, GCF_002859225.1_NZ_CP018761.1.region009, GCF_003186595.1_NZ_CP020645.1.region009, GCF_008693665.1_NZ_CP044089.1.region005, GCF_004006815.1_NZ_CP024175.1.region009, GCF_004008095.1_NZ_CP025069.1.region008
dim(distm):  2687 x 2687

=== Loading: syn.rds ===
str() summary:
List of 9
 $ clusterblast      : chr [1:21240] "APE Vf" "" "" "micacocidin" ...
 $ strainName        : chr [1:21240] "GCF_000009125.1" "GCF_000009125.1" "GCF_000009125.1" "GCF_000009125.1" ...
 $ regionName        : chr [1:21240] "NC_003295.1.region001" "NC_003295.1.region002" "NC_003295.1.region003" "NC_003295.1.region004" ...
 $ assemblyDefinition: chr [1:21240] "Ralstonia nicotianae GMI1000, complete sequence." "Ralstonia nicotianae GMI1000, complete sequence." "Ralstonia nicotianae GMI1000, complete sequence." "Ralstonia nicotianae GMI1000, complete sequence." ...
 $ location          : num [1:21240, 1:2] 433906 615572 1766202 1925813 3527740 ...
 $ wholeCDS          :List of 21240
 $ biosynCDS         :List of 21240
 $ group             : num [1:21240, 1] 0 0 0 51 0 0 0 33 0 0 ...
 $ regionIdentifier  : chr [1:21240] "GCF_000009125.1_NC_003295.1.region001" "GCF_000009125.1_NC_003295.1.region002" "GCF_000009125.1_NC_003295.1.region003" "GCF_000009125.1_NC_003295.1.region004" ...
NULL

Key head() checks:
head(strainName):  GCF_000009125.1, GCF_000009125.1, GCF_000009125.1, GCF_000009125.1, GCF_000009125.1, GCF_000009125.1
head(group):  0, 0, 0, 51, 0, 0

Inspection complete.

## Notes
- The `pairing_result.rds` object is a 26x2 matrix of lists, each containing a single element (likely the pairing data).
- `rec_cluster.rds` has a large distance matrix (42,190 x 42,190) for receptor clustering.
- `rec.rds` contains 42,190 receptor entries with detailed sequence and group information.
- `syn_cluster.rds` has a 2,687 x 2,687 distance matrix for synthetase clustering.
- `syn.rds` contains 21,240 synthetase entries, with many groups set to 0 (non-pyoverdine BGCs).
- All string fields were successfully decoded from HDF5 uint8 matrices.
- No loading errors occurred.
- For future inspections, re-run the script after any data updates.
