# -*- coding: utf-8 -*-


from coffe.core.placeholder import *


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
