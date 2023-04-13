## About Database
- The `gz` files were downloaded from NCBI ftp server on 2022-04-06.
- The data was bacteria 16S rRNA gene sequcnes (`fna`) and annotations (`gbff`).
- The `csv` file of type strinas -`tax-typestrain-20230406.csv`- was drawn from `gbff` file, containing taxonmy information:
	+ saccver = subject accestion version
	+ sname = subject name
	+ staxonomy = subject taxonomy
	+ sstrain = subjuect strain
	+ typestain = all True
- 1351 sequences that are not from type strain were deleted and saved as the clean fna file (`16S-type-20230406.fna`).
- The clean fna file was indexed by `makeblastdb`, and the output were 3 files: `16S.nhr`, `16S.nin` and `16s.nsp`.

## About data format
- The example of `gbff` is shown in `gbff-structure.gbff`.
- The example of `fna` is shown in `fna-example.fna`.
	+ When a fna file is imported by `SeqIO.parse()`, the fisrt word after `>` will be id and name, the whole line without `>` will be description.