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
    with open(str_file, 'w') as chmstr:
        chmstr.write(
            """* foo
*
read param card flex append

DIHEDRALS\n"""
        )
        for tfix in torfixes:
            atom1 = tfix.atom1
            atom2 = tfix.atom2
            atom3 = tfix.atom3
            atom4 = tfix.atom4
            value_dict = tfix.values
            for key, value in zip(value_dict.keys(), value_dict.values()):
                chmstr.write(
                    "{0:<7s}\t{1:<7s}\t{2:<7s}\t{3:<7s}\t{4:>7.4f}\t{5:>5d}\t{6:>6.1f}\n".format(
                        atom1.upper(), atom2.upper(),
                        atom3.upper(), atom4.upper(),
                        round((value[1]/4.184), 3),
                        int(key),
                        round((180*value[0]/3.141592653589793), 1)   
                    )
                )
        chmstr.write(
            """END
RETURN"""
        )

