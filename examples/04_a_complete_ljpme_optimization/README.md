# c36/ljpme

Workflow to optimize C36 lipid force field for LJPME simulation

indexing of properties:

0: area dppc bilayer @323.15K
1: area dppc bilayer @333.15K
2/3: DBs

Pre-requisite for new lipids/properties:

Define your lipids in $path_to_fflip/chm/lipids and add it to $path_to_fflip/chm/__init__.py
Create your template for the property calculation in $project_root/template (the root is core_files in this example)

After these:
(1) Run your simulation by: fflip simulate [OPTIONS]
(2) ...
