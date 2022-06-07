import argparse

parser = argparse.ArgumentParser(prog='gene_variants')
# parser.add_argument('--foo', action='store_true', help='foo help')
subparsers = parser.add_subparsers(help='sub-command help')

cmd_parser = subparsers.add_parser('parse', help='parse a VCF file')
cmd_parser.add_argument('file', type=str, help='path to the input VCF file')

args = parser.parse_args()

# # create the parser for the "b" command
# parser_b = subparsers.add_parser('b', help='b help')
# parser_b.add_argument('--baz', choices='XYZ', help='baz help')
