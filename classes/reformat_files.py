__author__ = 'davidmurphy'


from sys import argv
from shutil import copyfile, move


def reformat_fields(init_file):
    """Reformat init/final files RunStruct text files"""
    # copy file to backup before reformatting
    backup_file = init_file.replace('.txt', '.bkup.txt')
    copyfile(init_file, backup_file)

    # replace old fields with new fields throughout the file
    new_lines = []
    with open(init_file, 'r') as f:
        for line in f:

            # replacements
            line = line.replace('phase', 'phs')
            line = line.replace('token', 'tkn')
            line = line.replace('focal', 'fcl')
            line = line.replace('bkgd_scale', 'bscl')
            line = line.replace('bmap_dir', 'bdir')
            line = line.replace('cmap_dir', 'cdir')
            line = line.replace('bs_dfe', 'bdfe')
            line = line.replace('cs_dfe', 'cdfe')
            line = line.replace('label', 'lbl')
            # removals
            if 'rm_annos' in line:
                continue

            new_lines.append(line)

    with open(init_file, 'w') as f:
        f.write(''.join(line for line in new_lines))

    return None


def rename_files(file_name, old_pattern, new_pattern):
    """Rename file by replacing an old pattern with a new pattern"""
    # ntarr -> narr
    # neut/mask -> nmsk
    # neut/mask/rand -> nmsk/rdm


def main():
    if len(argv) != 2:
        print 'usage: reformat_files <init_file>'
        exit(1)

    init_file = argv[1]
    reformat_fields(init_file)


if __name__ == '__main__':
    main()
