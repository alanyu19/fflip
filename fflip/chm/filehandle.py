# -*- coding: utf-8 -*-


from coffe.core.placeholder import *
from coffe.omm.torsionfuncs import *


def replace_top_parameter(top_file_template, para_dict, write_to_new_dir=None):
    if write_to_new_dir is None:
        replace_placeholders(top_file_template, para_dict=para_dict)
    else:
        assert os.path.isdir(write_to_new_dir)
        os.system('cp {} {}'.format(top_file_template, write_to_new_dir))
        new_file = os.path.join(write_to_new_dir, top_file_template.split('/')[-1])
        replace_top_parameter(
            new_file,
            para_dict=para_dict
        )


def torfix_to_str(torfixes, str_file):
    with open(str_file, 'w') as str:
        for tfix in torfixes:
            atom1 = tfix.atom1
            atom2 = tfix.atom2
            atom3 = tfix.atom3
            atom4 = tfix.atom4
            value_dict = tfix.values
            for key, value in zip(value_dict.keys, value_dict.values):
                str.write(
                    atom1.upper(), atom2.upper(), atom3.upper(), atom4.upper(),
                    str(key), str(value[0]), str(value[1])
                )

