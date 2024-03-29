NGS Scripts
===========

**Jason Ross <jason.ross at csiro.au>**

compile_fastqc_data.py
----------------------


Classes for parsing FASTQC reports into Pandas objects and for collating a series of FASTQC reports in a directory tree into one report per FASTQC module.
The *fastqc_collation* class will collate all FASTQC reports within a folder branch specified by project_root into a report for each FASTQC module. The output directory, FASTQC filename and field seperator in the report can optionally be specified.
The *fastqc_report* class holds an instance of one FASTQC report. The modules are contained in a dict of *fastqc_module*, which itself is a basic Pandas container class.

For basic usage, the module may be run from the command line by specifying a project root folder and optionally an output directory (defaults to project_root):

	compile_fastqc_data.py --project_root '/my/project/folder' --out_dir '/my/out/folder'

The classes may also be imported and used like below:

	>>>project_root = '/my/project/folder'
	>>> collation = fastqc_collation(project_root, out_dir = '/my/out/folder')
	>>> collation.modules
	['Per base sequence quality', 'Sequence Duplication Levels', 'Kmer Content', 'Per base sequence content', 'Per sequence GC content', 'Sequence Length Distribution', 'Per base GC content', 'Basic Statistics', 'Overrepresented sequences', 'Per base N content', 'Per sequence quality scores']

For for just the 'Per sequence GC content' module output, we type:

	collation.process('Per sequence GC content')

For all modules to be processed, type:

	collation.process()

The fastqc_report class can be used to process a single FASTQC file:

	>>> fastqc = fastqc_report('/my/project/folder/subfolder/fastqc_data.txt')
	>>> fastqc['Kmer Content']
	   Sequence    Count Obs/Exp Overall Obs/Exp Max Max Obs/Exp Position
	0     TGAAA  4567880        2.992366     6.29317                    2
	1     AAAAA  4099530       2.8870537    8.209245                60-64
	2     TCAAA  4183840        2.726956   5.4784513                    7
	3     TTTGA  3981485       2.5597813   5.1077914                    6
	4     AAATT  3574600       2.4706206   5.9077806                    4
	5     GTGAA  3767510         2.31742    5.057868                    1
	(and so on)




