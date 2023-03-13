import datetime
import gzip
import hashlib
import multiprocessing
import os
import re

import config
from database import Database


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
                                          args=(file_path, file_chunk_size, db_chunk_size), )
                results.append(result)

            for result in results:
                result.wait()

            end_time = datetime.datetime.now().replace(microsecond=0)
            print(f'Time taken to process all files = {end_time - start_time}')
        else:
            raise RuntimeError(f"Error: invalid file type '{file_type}'")


def _compress_compute_hash(file_path, file_chunk_size):
    hasher = hashlib.md5()
    file_name = os.path.basename(file_path)
    temp_file_path = os.path.join(config.FASTA_FILE_OUTPUT_DIRECTORY, f"{file_name}.gz")

    if not os.path.exists(os.path.dirname(temp_file_path)):
        os.makedirs(os.path.dirname(temp_file_path))

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
    # share an object between processes using multiprocessing.Manager().Namespace() (check pickling)
    db = Database()

    if not os.path.isfile(file_path):
        raise ValueError(f'File not found: {file_path}')

    file_name = os.path.basename(file_path)
    start_time = datetime.datetime.now().replace(microsecond=0)

    print(f'Processing file {file_name} on {multiprocessing.current_process().name}')
    final_file_path, file_checksum = _compress_compute_hash(file_path=file_path, file_chunk_size=file_chunk_size)

    print(f'Loading file {os.path.basename(final_file_path)} to db')

    # check if file is already processed
    if _find_file_by_name_checksum(db=db, file_name=file_name, file_checksum=file_checksum):
        print(f'File is already processed. file_name={file_name}, file_checksum={file_checksum}')
        return

    file_type_format = 'text/x-fasta'
    file_doc = _create_file_document(file_name=file_name, file_checksum=file_checksum,
                                     db=db, file_type_format=file_type_format)
    print(f"file_doc inserted = {file_doc.inserted_id}")

    next_sequence_val = _get_next_sequence_value(db=db)
    assembly_name = os.path.splitext(os.path.basename(file_name))[0]
    user = 'python-loader-' + str(next_sequence_val['sequenceValue'])

    # check if assembly exists in db
    if _find_assembly_by_name(db=db, assembly_name=assembly_name):
        print(f'Assembly already exists. assembly_name={assembly_name}')
        return

    # create assembly document
    assembly_doc = _create_assembly_document(db=db, assembly_name=assembly_name, user=user)
    assembly_id = assembly_doc.inserted_id
    print(f"assembly_doc inserted = {assembly_id}")

    # read the sequences and write to db
    try:
        sequence = ''
        chunk_index = 0
        ref_seq_len = 0
        ref_seq_doc = None

        with gzip.open(final_file_path, 'rb') as fi:
            for line in fi:
                decoded_line = line.rstrip().decode()

                # extract sequence name and description from header
                ref_seq_info_line = re.match(r'^>\s*(\S+)\s*(.*)', decoded_line)

                # parse if its header line
                if ref_seq_info_line:
                    # If there is unsaved sequence of previous sequence header
                    if sequence != '':
                        # if header line not present
                        if not ref_seq_doc:
                            raise Exception('No refSeq document found')

                        ref_seq_len += len(sequence)
                        print(f'refseqchunks inserted for refseqs id = {ref_seq_doc.inserted_id}')
                        _create_refseqchunks(db=db, inserted_id=ref_seq_doc.inserted_id, index=chunk_index,
                                             sequence=sequence, user=user)
                        sequence = ''

                    if ref_seq_doc:
                        # update the ref sequence length
                        _update_refseqs_len(db=db, inserted_id=ref_seq_doc.inserted_id, asm_id=assembly_id,
                                            length=ref_seq_len)

                    # Re-initialize chunk index and ref_seq_len to 0 as we are processing different sequence header
                    chunk_index = 0
                    ref_seq_len = 0
                    # get sequence name and description from regex groups
                    ref_seq_name = ref_seq_info_line.group(1)
                    ref_seq_desc = ref_seq_info_line.group(2)
                    ref_seq_doc = _create_refseqs(db=db, name=ref_seq_name, desc=ref_seq_desc, asm_id=assembly_id,
                                                  db_chunk_size=db_chunk_size, user=user)
                    print(f'refseqs inserted = {ref_seq_doc.inserted_id}')
                else:
                    # if header line not present
                    if not ref_seq_doc:
                        raise Exception('No refSeq document found')

                    if re.search(r'\S', decoded_line):
                        sequence += decoded_line
                        while len(sequence) >= db_chunk_size:
                            sequence_chunk = sequence[:db_chunk_size]
                            ref_seq_len += len(sequence_chunk)
                            print(f'refseqchunks inserted for refseqs id = {ref_seq_doc.inserted_id}')
                            _create_refseqchunks(db=db, inserted_id=ref_seq_doc.inserted_id, index=chunk_index,
                                                 sequence=sequence_chunk, user=user)
                            chunk_index += 1
                            sequence = sequence[db_chunk_size:]

        # if there are remaining sequence to be saved, happens when sequence size is < db_chunk_size
        if sequence != '':
            ref_seq_len += len(sequence)
            _create_refseqchunks(db=db, inserted_id=ref_seq_doc.inserted_id, index=chunk_index,
                                 sequence=sequence, user=user)
            print(f'remaining refseqchunks inserted for refseqs id = {ref_seq_doc.inserted_id}')
            _update_refseqs_len(db=db, inserted_id=ref_seq_doc.inserted_id, asm_id=assembly_id,
                                length=ref_seq_len)

    except Exception as ex:
        print(ex)

    change_doc = _create_change_document(assembly_id=assembly_id, file_name=file_name, db=db, file_doc=file_doc,
                                         next_sequence_val=next_sequence_val)
    print(f"change_doc inserted = {change_doc.inserted_id}")

    end_time = datetime.datetime.now().replace(microsecond=0)
    print(f'Time taken to process file {file_name} = {end_time - start_time}')


def _create_refseqs(db, name, desc, asm_id, db_chunk_size, user):
    return db.refseqs.insert_one({
        'name': name,
        'description': desc,
        'assembly': asm_id,
        'length': 0,
        'chunkSize': db_chunk_size,
        'user': user,
        'status': 0,
    })


def _update_refseqs_len(db, inserted_id, asm_id, length):
    db.refseqs.update_one({'_id': inserted_id, 'assembly': asm_id},
                          {'$set': {'length': length}})


def _create_refseqchunks(db, inserted_id, index, sequence, user):
    db.refseqchunks.insert_one({
        'refSeq': inserted_id,
        'n': index,
        'sequence': sequence,
        'user': user,
        'status': 0,
    })


def _get_next_sequence_value(db):
    return db.counters.find_one_and_update(
        {'id': 'changeCounter'},
        {'$inc': {'sequenceValue': 1}},
        upsert=True,
        return_document=True
    )


def _create_change_document(assembly_id, file_name, db, file_doc, next_sequence_val):
    change = {
        'typeName': 'AddAssemblyFromFileChange',
        'assembly': str(assembly_id),
        'assemblyName': os.path.splitext(os.path.basename(file_name))[0],
        'fileId': str(file_doc.inserted_id)
    }
    change_doc = {
        'typeName': 'AddAssemblyFromFileChange',
        'assembly': assembly_id,
        'changes': [change],
        'user': 'python-loader-' + str(next_sequence_val['sequenceValue']),
        'sequence': next_sequence_val['sequenceValue']
    }
    return db.changes.insert_one(change_doc)


def _find_file_by_name_checksum(db, file_name, file_checksum):
    try:
        return db.files.find_one({
            'basename': file_name,
            'checksum': file_checksum,
        })
    except Exception as ex:
        print(f'Error at _find_file_by_checksum. {ex}')
        return None


def _create_file_document(file_name, file_checksum, db, file_type_format):
    return db.files.insert_one({
        'basename': file_name,
        'checksum': file_checksum,
        'type': file_type_format
    })


def _create_assembly_document(db, assembly_name, user):
    return db.assemblies.insert_one({
        'name': assembly_name,
        'user': user,
        'status': 0
    })


def _find_assembly_by_name(db, assembly_name):
    try:
        return db.assemblies.find_one({
            'name': assembly_name
        })
    except Exception as ex:
        print(f'Error at _find_assembly_by_name. {ex}')
        return None
