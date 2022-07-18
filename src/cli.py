import argparse
import sys

import tasks

parser = argparse.ArgumentParser(prog='gene_variants')
# parser.add_argument('--foo', action='store_true', help='foo help')
subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')

cmd_parser = subparsers.add_parser('parse', help='parse a VCF file')
cmd_parser.add_argument('vcf_file', type=str, help='path to the input VCF file')
cmd_parser.add_argument('genes_file', type=str, help='path to a file that includes one gene of interest per line (as an HGNC)')

# # create the parser for the "b" command
# parser_b = subparsers.add_parser('b', help='b help')
# parser_b.add_argument('--baz', choices='XYZ', help='baz help')

# print help when no arguments are provided
if len(sys.argv) == 1:
    parser.print_help()

args = parser.parse_args()

if args.subcommand == 'parse':
    # TODO: save and pass gene set properly
    tasks.parse(args.vcf_file, args.genes_file, -1)
else:
    print('no can do')
    exit(1)