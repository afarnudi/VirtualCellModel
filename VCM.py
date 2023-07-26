import sys

from source.general.argument_parser import analyse_parser
from source.general.argument_parser import create_parser

def main():
    parser = create_parser()
    user_args = parser.parse_args()
    user_inputs = analyse_parser(user_args, parser)


if __name__ == "__main__":
    main()
