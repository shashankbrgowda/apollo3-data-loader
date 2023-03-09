import datetime
import gzip
import hashlib
import multiprocessing
import os

import config
from database import Database

COMPUTE_AND_COMPRESS_FILE = True


def _compress_compute_hash(file_path, file_chunk_size):
    hasher = hashlib.md5()
    file_name = os.path.basename(file_path)
    temp_file_path = os.path.join(config.FASTA_FILE_OUTPUT_DIRECTORY, f"{file_name}.gz")

    try:
        with open(file_path, 'rb') as fi, gzip.open(temp_file_path, 'wb') as fo:
            # Read data from the file in chunks of a specified chunk size
            for data_chunk in iter(lambda: fi.read(file_chunk_size), b''):
                # Update the hash object with each chunk of data
                hasher.update(data_chunk)
                fo.write(data_chunk)

            # Compute the hash digest of the entire file
            file_checksum = hasher.hexdigest()
            final_file_path = os.path.join(config.FASTA_FILE_OUTPUT_DIRECTORY, file_checksum)
            os.rename(temp_file_path, final_file_path)
            print(f'File name = {file_name}, checksum = {file_checksum}, compressed file path = {final_file_path}')
            return final_file_path, file_checksum
    except Exception as ex:
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)
        print(f'Error: failed to read fasta file = {file_name}, error = {ex}')


def _process_fasta_file(file_path, file_chunk_size, db_chunk_size):
    """
    Compresses a file using gzip, saves it to a file with a name based on its checksum and reads the
    sequence from the compressed final to write it to mongodb
    """

    # connection object is not thread safe, it's better to have a connection object per process
    # share an object between processes using multiprocessing.Manager().Namespace() but, its tricky (Pickling)
    db = Database()

    if not os.path.isfile(file_path):
        raise ValueError(f'File not found: {file_path}')

    file_name = os.path.basename(file_path)
    start_time = datetime.datetime.now().replace(microsecond=0)

    print(f'Processing file {file_name} on {multiprocessing.current_process().name}')
    final_file_path, file_checksum = _compress_compute_hash(file_path=file_path, file_chunk_size=file_chunk_size)

    if COMPUTE_AND_COMPRESS_FILE:
        print(f'Loading file {os.path.basename(final_file_path)} to db')
        print(db.assemblies.find_one())
    else:
        print(f'Loading file {file_name} to db')

    end_time = datetime.datetime.now().replace(microsecond=0)
    print(f'Time taken to process file {file_name} = {end_time - start_time}')


def process(file_chunk_size, db_chunk_size, file_type):
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        if file_type == config.FASTA_FILE_TYPE:
            fasta_files = os.listdir(config.FASTA_FILES_DIRECTORY)
            print(f'Fasta files to be processed = {fasta_files}')

            results = []
            start_time = datetime.datetime.now().replace(microsecond=0)
            for fasta_file in fasta_files:
                file_path = os.path.join(config.FASTA_FILES_DIRECTORY, fasta_file)
                result = pool.apply_async(_process_fasta_file,
                                          args=(file_path, file_chunk_size, db_chunk_size),)
                results.append(result)

            for result in results:
                result.wait()

            end_time = datetime.datetime.now().replace(microsecond=0)
            print(f'Time taken to process all files = {end_time - start_time}')
        else:
            raise RuntimeError(f'Error: invalid file type \'{file_type}\'')
