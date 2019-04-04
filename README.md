## VCF indel notation converter from (A -> "-") to (GA -> G) including leading base (proper VCF notation).

#### For info please contact me: or.yaacov@mail.huji.ac.il

##### Dependencies:
R (3+), 
Bioconductor (packages: BSgenomem, BSgenome.Hsapiens.UCSC.hg19, Biostrings)

Takes a 5 col tsv file (chr, pos, name, ref, alt):
>chr1	20996757	NULL	T	- <br>
>chr1	20996257	NULL	TT	- <br>
>chr1	20996457	NULL	-	TT <br>
>chr1	20996457	NULL	-	T <br>
>chr1	20996457	NULL	A	G<br>

Converts to:
>chr1	20996756	NULL	AT	A<br>
>chr1	20996255	NULL	AGTT	AG<br>
>chr1	20996455	NULL	GA	GATT<br>
>chr1	20996456	NULL	A	AT<br>
>chr1	20996457	NULL	A	G<br>
