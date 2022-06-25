import hashlib
import os


def each_cons(xs, n):
    return [xs[i:i+n] for i in range(len(xs)-n+1)]


def sha256sum(filename):
    h  = hashlib.sha256()
    b  = bytearray(128*1024)
    mv = memoryview(b)
    with open(filename, 'rb', buffering=0) as f:
        for n in iter(lambda : f.readinto(mv), 0):
            h.update(mv[:n])
    return h.hexdigest()


def get_data_dir(file):
    dir_name = os.path.basename(file)
    return 'data/intermediary/{}'.format(dir_name)


def get_genes_from_file(file):
    return [gene.decode('utf-8').rstrip() for gene in file]