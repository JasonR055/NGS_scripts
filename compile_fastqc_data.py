#!/usr/bin/env python
"""
    compile_fastqc_data
    Jason Ross <jason.ross@csiro.au>
    
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
    
"""

import os, sys, re, argparse, errno
import pandas as pd
from numpy import array


class fastqc_collation(object):
    """Collate all FASTQC reports within a folder branch specified by
    project_root into a report for each FASTQC module. The output directory,
    FASTQC filename and field seperator in the report can optionally be
    specified.
    """
    
    def __init__(self, project_root, out_dir = None,
                 fastqc_filename='fastqc_data.txt', sep=','):
        
        if out_dir is None:
            out_dir = project_root
        self.project_root = project_root
        self.out_dir = out_dir
        self.fastqc_filename = fastqc_filename
        self.sep = sep
        self.fastqcfiles = []
        
        for root, dirs, files in os.walk(self.project_root):
            # walk a directory containing FastQC output for multiple samples
            for name in files:
                if (name == fastqc_filename):
                    fastqc_fullpath = os.path.join(root, fastqc_filename)
                    self.fastqcfiles.append(fastqc_fullpath)
        self.fastqcfiles.sort()
        # Get the exhaustive list of modules from the first report
        self.modules = fastqc_report(self.fastqcfiles[0]).modules.keys()
    
    
    def process(self, modules=None):
        """Process the given modules in the collated FASTQC files in turn. If
        no modules are specified (default), all are processed).
        
        Example:
            process('Per sequence GC content')
        """
        
        out_files = {}
        self._make_sure_path_exists()
        suffix = '.csv'
        if self.sep != ',':
            suffix='.txt'
        
        if modules is None:  # Process all modules
            modules = self.modules
        
        if type(modules) is str:
            modules = [modules]
        
        for module in modules:
            out_name = module + '_fastqc_collation' + suffix
            out_files[module] = os.path.join(self.out_dir, out_name)
            # Remove old files if present
            if os.path.isfile(out_files[module]):
                os.remove(out_files[module])
        
        for f in self.fastqcfiles:
            qc = fastqc_report(f)
            for module in modules:
                # If the file exists, append it, otherwise create it
                if(os.path.isfile(out_files[module])):
                    self._process_mod(f, qc, module, out_files[module],
                                      mode='a')
                else:
                    self._process_mod(f, qc, module, out_files[module])
        return
    
    
    def _process_mod(self, filepath, qc_report, module, out, mode='w'):
        """Write or append to a module file"""
        
        if module not in self.modules:
            raise KeyError("Module not found.")
        if mode not in ['w', 'a']:
            raise ValueError("Only write (w) and append (a) modes accepted.")
        
        data = qc_report.modules[module].data
        data['Filestub'] = qc_report.filestub
        data['Filepath'] = filepath
        
        if module == 'Basic Statistics':
            fh = open(out, mode=mode)
            if mode == 'w':
                # Create a header line
                fh.write(self.sep.join(list(data.index)))
                fh.write('\n')
            fh.write(self.sep.join(data))
            fh.write('\n')
            fh.close()
        else:
            if mode == 'w':
                header=True
            else:
                header=False
            data.to_csv(out, index=False, sep=self.sep,
                        header=header, mode=mode)
        return
        
    
    def _make_sure_path_exists(self):
        try:
            os.makedirs(self.out_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


class fastqc_report(object):
    """Process a FASTQC report file. Each module is contained by a
    fastqc_module object."""
    
    def __init__(self, fastqcfile):
        
        self._check_file_exists(fastqcfile)
        self.fastqcpath = fastqcfile
        self.modules = {}
        self._process_fastqc_report()
        self._process_basic_stats()
    
    def _check_file_exists(self, fastqc_path):
        if not os.path.isfile(fastqc_path):
            raise ValueError("File doesn't exist")
    
    def _process_fastqc_report(self):
        
        lines = []
        with open(self.fastqcpath, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if (line[:2] == ">>" and line[:12] != ">>END_MODULE"):
                    [module, status] = line[2:].split('\t')
                    continue
                if (line[:1] == "#"):
                    cols = line[1:].split('\t')
                    continue
                if line == ">>END_MODULE":
                    mod = fastqc_module(module, status, cols, lines)
                    self.modules[module] = mod
                    lines = []
                else:
                    lines.append(line.split('\t'))
    
    def _process_basic_stats(self):
        stats = self.modules['Basic Statistics'].data
        self.seq_total = int(stats['Total Sequences'])
        self.seq_filtered = int(stats['Filtered Sequences'])
        self.seq_length = int(stats['Sequence length'])
        self.seq_gc = int(stats['%GC'])
        self.filename = stats['Filename']
        self.filestub = re.sub('(.*)\.fastq.*', '\\1', self.filename)

    
    def __getitem__(self, key):
        return self.modules[key].data


class fastqc_module(object):
    """Holds FASTQC module data in a Pandas object"""
    
    def __init__(self, module, status, cols, lines):
        
        self.module = module
        self.status = status
        self.cols = cols
        
        if module == 'Basic Statistics':
            lines_array = array(lines)
            self.data = pd.Series(lines_array[:,1], index=lines_array[:,0])
        else:
            self.data = pd.DataFrame(lines, columns=cols)
    
    #def __repr__(self):
    #    return str(self.data)


if __name__ == '__main__':
	
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--project_root', type=str, default=sys.stdin)
    parser.add_argument('-o', '--out_dir', nargs='?', type=str, default=None)
    parser.add_argument('-f', '--fastqc_filename', nargs='?', type=str,
                        default='fastqc_data.txt')
    parser.add_argument('-s', '--sep', nargs='?', type=str, default=',')
    ns = parser.parse_args()
    collation = fastqc_collation(ns.project_root, ns.out_dir,
                                 ns.fastqc_filename, ns.sep)
    collation.process()
