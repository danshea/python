#!/usr/bin/env python
#
# Name: Galaxy2Cummerbund.py
# Date: 2014-06-30
# Author: Dan Shea <daniel.john.shea@gmail.com>
# Description:
#
# Galaxy does not use the same file naming convention as cummeRbund
# So, we need to rename the files
#
# Galaxy Name                                       cummeRbund Name
#----------------------------------------------------------------------------
# transcript_FPKM_tracking                          isoforms.fpkm_tracking
# transcript_differential_expression_testing        isoform_exp.diff
# gene_FPKM_tracking                                genes.fpkm_tracking
# gene_differential_expression_testing              gene_exp.diff
# TSS_groups_FPKM_tracking                          tss_groups.fpkm_tracking
# TSS_groups_differential_expression_testing        tss_group_exp.diff
# CDS_FPKM_tracking                                 cds.fpkm_tracking
# CDS_FPKM_differential_expression_testing          cds_exp.diff
# CDS_overloading_differential_expression_testing   cds.diff
# promoters_differential_expression_testing         promoters.diff
# splicing_differential_expression_testing          splicing.diff

import argparse
import os
import os.path
import re
import shutil
import sys

def rename_files(directory):
    conversion = {
                'transcript_FPKM_tracking': 'isoforms.fpkm_tracking',
                'transcript_differential_expression_testing': 'isoform_exp.diff',
                'gene_FPKM_tracking': 'genes.fpkm_tracking',
                'gene_differential_expression_testing': 'gene_exp.diff',
                'TSS_groups_FPKM_tracking': 'tss_groups.fpkm_tracking',
                'TSS_groups_differential_expression_testing': 'tss_group_exp.diff',
                'CDS_FPKM_tracking': 'cds.fpkm_tracking',
                'CDS_FPKM_differential_expression_testing': 'cds_exp.diff',
                'CDS_overloading_differential_expression_testing': 'cds.diff',
                'promoters_differential_expression_testing': 'promoters.diff',
                'splicing_differential_expression_testing': 'splicing.diff',
    }
    for filename in os.listdir(directory):
        for key in conversion.keys():
            pattern = '.*({})'.format(key)
            match = re.match(pattern, filename)
            if match is not None:
                lookup = match.group(1)
                dst = os.path.join(directory,conversion[lookup])
                src = os.path.join(directory,filename)
                print 'Renaming {} to {}'.format(src, dst)
                shutil.copy2(src, dst)

def main():
    parser = argparse.ArgumentParser(description='Covert Galaxy cuffdiff files to the expected names for cummeRbund analysis.')
    parser.add_argument('directory',help='directory where Galaxy files are located')
    args = parser.parse_args()
    print 'directory is {}'.format(args.directory)
    rename_files(args.directory)
    return(0)

if __name__ == '__main__':
    exit(main())
