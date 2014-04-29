NGS Scripts
===========

**Jason Ross <jason.ross at csiro.au>**

compile_fastqc_data
-------------------


Classes for parsing FASTQC reports into Pandas objects and for collating
a series of FASTQC reports in a directory tree into one report per FASTQC
module.

The module may be run from the command line or imported and used like
below::

	project_root = '/my/project/folder'
	collation = fastqc_collation(project_root, out_dir = '/my/out/folder')

For for just the 'Per sequence GC content' module output, we type::

	collation.process('Per sequence GC content')

For all modules to be processed, type::

	collation.process()



