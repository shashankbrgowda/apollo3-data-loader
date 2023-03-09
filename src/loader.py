import sys

import processor
import config

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def _argument_parser() -> ArgumentParser:
    parser = ArgumentParser(description='Apollo3 mongodb data loader',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--file-type', type=str, required=True, help='fasta or gff')
    parser.add_argument('--file-chunk-size', type=int, help='file chunk size')
    parser.add_argument('--db-chunk-size', type=int, help='db chunk size')
    return parser


if __name__ == '__main__':
    parser = _argument_parser()
    args = parser.parse_args()
    print(f'Passed arguments = {vars(args)}')

    file_type = args.file_type
    file_chunk_size = args.file_chunk_size if args.file_chunk_size is not None else config.FILE_CHUNK_SIZE
    db_chunk_size = args.db_chunk_size if args.db_chunk_size is not None else config.DB_CHUNK_SIZE

    try:
        processor.process(file_chunk_size=file_chunk_size,
                          db_chunk_size=db_chunk_size,
                          file_type=file_type)
    except RuntimeError as ex:
        print(ex)
        sys.exit(1)
