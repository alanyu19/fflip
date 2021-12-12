# -*- coding: utf-8 -*-

import math


class nbgroup(object):
    """
    A class that contains the info about the group for which one want to
    change the nobonded parameters but not the interaction form ...
    "nbgroup" = GroupToChangeNonbondedParameterss
    """
    def __init__(self, **kwargs):
        self.par_type = kwargs["par_type"]
        self.center_names = kwargs["center_names"]
        if self.par_type == "charge":
            assert "cooperators" in kwargs
            self.cooperators = kwargs["cooperators"]
            assert "neighbors" in kwargs
            self.neighbors = kwargs["neighbors"]

        # keep in mind that the "original_p" should
        # be the same as those in the parameter file
        # (at least for the current version of this code)
        self.original_p = kwargs["original_p"]
        # default offset might be useful in energy evaluation...
        if "default_offset" in kwargs:
            self.default_offset = kwargs["default_offset"]

        # there are two ways to play around with the parameter in optimization,
        # either give a list of targeted parameters/scale_factors
        # OR
        # give the offset range directly
        # we are currently using the first scheme since
        # it's more straightforward ...
        # self.offset_range = kwargs["offset_range"]
        # note that targeted ranges in this version of code are:
        # 1. REAL parameter ranges for charges;
        # 2. range of scale_factors for epsilon and sigma.
        if self.par_type == "charge" or (self.par_type != "charge"):
            self.targeted_range = kwargs["targeted_range"]
            # TODO: convert to scale functionality (YYL)

        if self.par_type == "charge":
            # "roc" stands for the "ratios of cooperators
            # (to center in terms of change of charges)"
            assert "roc" in kwargs
            self.roc = kwargs["roc"]
            assert "ron" in kwargs
            self.ron = kwargs["ron"]
            # We might want to add a checker here to verify
            # the net charge of the group is conserved.

    @property
    def offset_range(self):
        """
        find the offset for the global parameters, might need some change
        in the pre-factors depending on the units
        :return:
        """
        offset = []
        if self.par_type == "charge":
            offset.append(round(self.targeted_range[0] - self.original_p, 4))
            offset.append(round(self.targeted_range[1] - self.original_p, 4))
        elif self.par_type == "sigma" or self.par_type == "epsilon":
            # epsilon or sigma
            # this might not be useful if the default offset is
            # exactly the original parameter
            offset.append(
                round((self.targeted_range[0] - 1) * self.original_p, 4)
            )
            offset.append(
                round((self.targeted_range[1] - 1) * self.original_p, 4)
            )
        else:
            offset = [0, 0]  # is pass a better option?
        return offset

    @property
    def scale_offset_range(self):
        """
        this is useful if the per-particle offset for sigma and epsilon
        can be set to be the original parameter. In the worst case, we can
        keep using offset_range above.
        """
        scale_offset = []
        if self.par_type == "sigma" or self.par_type == "epsilon":
            scale_offset.append(round(self.targeted_range[0] - 1, 2))
            scale_offset.append(round(self.targeted_range[1] - 1, 2))
        # in most cases, this will be a negative value followed by
        # a positive value.
        return scale_offset

    @property
    def force_names(self):
        name = self.center_names[0] + '_' + self.par_type
        if self.par_type == "sigma" or self.par_type == "epsilon":
            return [name]
        elif self.par_type == "charge":
            # initialize name list
            name_list = [name]
            coop_names = []
            for coop in self.cooperators:
                coop_names.append(
                    coop + '_' + self.center_names[0] + '_' + self.par_type
                )
                name_list.append(
                    coop + '_' + self.center_names[0] + '_' + self.par_type
                )
            neib_names = []
            for neib in self.neighbors:
                neib_names.append(
                    neib + '_' + self.center_names[0] + '_' + self.par_type
                )
                name_list.append(
                    neib + '_' + self.center_names[0] + '_' + self.par_type
                )
            return [name, coop_names, neib_names]
        else:
            return None

    def __repr__(self):
        # par_type, center_names,
        if self.par_type == "charge":
            return "nbgroup(par_type = '{}', center_names = {}," \
                "original_p = {}, targeted_range = {}, neighbors = {}," \
                " cooperators = {}, ron = {}, roc = {})\n".format(
                    self.par_type, self.center_names, self.original_p,
                    self.targeted_range, self.neighbors, self.cooperators,
                    self.ron, self.roc
                )
        elif self.par_type == "epsilon" or self.par_type == "sigma":
            return "nbgroup(par_type = '{}', center_names = {}," \
                "original_p = {}, targeted_range = {})\n".format(
                    self.par_type, self.center_names,
                    self.original_p, self.targeted_range
                )
        else:
            return "nbgroup(no change)\n"


class GTCTS:
    """
    a class that contains the info about the group for which one want to
    change the charmm torsion force constant and potentially multiplicity ...
    "gtctp" = GroupToChangeTorsionParameters
    """
    # phase two of my work (maybe not needed at all ...)


def empty_nbgroup():
    return nbgroup(
        par_type=None, center_names=[], original_p=0, targeted_range=[]
    )

# nbgroup(par_type="charge", center_names=["N"], cooperators=[],
# neighbors=['C12', 'C13', 'C14', 'C15'], original_p=0.5, targeted_range=[
# 0.4, 0.7], roc=[], ron=[-0.25, -0.25, -0.25, -0.25])

# -----------------------------------------------------------------------------
# DRUDE -----------------------------------------------------------------------


class DrudeParameter:
    def __init__(self, lipidname, cgid, internal_id,
                 par_type, center_names, original_p, targeted_range,
                 neighbors=None, drude_particles=None, **kwargs):
        """
        Class contains
        Args:
            lipidname: str, the lipid it belongs to
            cgid: int, the group id it belongs to (see dlipid)
            partype: can be 'charge, alpha, thole, sigma, epsilon
            center_names:
            neighbors:
            original_p:
            targeted_range:
            **kwargs:
        """
        self.lipid = lipidname
        self.cgid = cgid
        self.internal_id = internal_id
        self.par_type = par_type
        self.center_names = center_names
        if self.par_type == "charge":
            assert neighbors is not None
            self.neighbors = neighbors
        if self.par_type == "alpha" or self.par_type == "thole":
            assert drude_particles is not None  # maybe not needed
            self.drude_particles = drude_particles

        # keep in mind that the "original_p" should
        # be the same as those in the parameter file
        # (at least for the current version of this code)
        self.original_p = original_p
        # default offset might be useful in energy evaluation...
        if "default_offset" in kwargs:
            self.default_offset = kwargs["default_offset"]

        # there are two ways to play around with the parameter in optimization,
        # either give a list of targeted parameters/scale_factors
        # OR
        # give the offset range directly
        # we are currently using the first scheme since
        # it's more straightforward ...
        # self.offset_range = kwargs["offset_range"]
        # note that targeted ranges in this version of code are:
        # 1. REAL parameter ranges for charges;
        # 2. range of scale_factors for epsilon and sigma.

        self.targeted_range = targeted_range  # how did I use this ?????

        if self.par_type == "charge":
            assert "ron" in kwargs
            self.ron = kwargs["ron"]

    @property
    def offset_range(self):
        """
        Offset for the global parameters, might need some change
        in the pre-factors depending on the units.
        """
        offset = []
        if self.par_type == "charge" or self.par_type == 'alpha' \
            or self.par_type == 'thole' or self.par_type == 'nbthole':
            offset.append(
                round(self.targeted_range[0] - self.original_p, 4)
            )
            offset.append(
                round(self.targeted_range[1] - self.original_p, 4)
            )
        elif self.par_type == "sigma" or self.par_type == "epsilon":
            # epsilon or sigma
            # this might not be useful if the default offset is
            # exactly the original parameter
            offset.append(
                round((self.targeted_range[0] - 1) * self.original_p, 4)
            )
            offset.append(
                round((self.targeted_range[1] - 1) * self.original_p, 4)
            )
        else:
            offset = [0, 0]  # no change, for original energy calculation
        return offset

    @property
    def scale_offset_range(self):
        """
        This is useful if the per-particle offset for sigma and epsilon
        can be set to be the original parameter. In the worst case, we can
        keep using offset_range above.
        """
        scale_offset = []
        if self.par_type == "sigma" or self.par_type == "epsilon":
            scale_offset.append(round(self.targeted_range[0] - 1, 2))
            scale_offset.append(round(self.targeted_range[1] - 1, 2))
        # in most cases, this will be a negative value followed by
        # a positive value.
        return scale_offset

    @property
    def force_names(self):
        # this is not useful anymore probably
        name = self.center_names[0] + '_' + self.par_type
        if self.par_type == "sigma" or self.par_type == "epsilon" \
                or self.par_type == "alpha" or self.par_type == "thole" \
                or self.par_type == "nbthole":
            return [name]
        elif self.par_type == "charge":
            # initialize name list
            neighbor_names = []
            for ngb in self.neighbors:
                neighbor_names.append(
                    ngb + '_' + self.center_names[0] + '_' + self.par_type
                )
            return [name, neighbor_names]
        else:
            return None

    def __repr__(self):
        # par_type, center_names,
        if self.par_type == "charge":
            return "DrudeParameter of {} group {} (par_type='{}', " \
                   "center_names={}, original_p={}, targeted_range={}, " \
                   "neighbors={}, ron={})\n".format(
                        self.lipid, self.cgid,
                        self.par_type, self.center_names, self.original_p,
                        self.targeted_range, self.neighbors, self.ron
                    )
        elif self.par_type == "alpha" or self.par_type == 'thole':
            return "DrudeParameter of {} group {} (par_type='{}', " \
                   "center_names={}, drude_particles={}, original_p={}, " \
                   "targeted_range={}\n".format(
                        self.lipid, self.cgid,
                        self.par_type, self.center_names, self.drude_particles,
                        self.original_p, self.targeted_range
                    )
        elif self.par_type == "epsilon" or self.par_type == "sigma":
            return "DrudeParameter(par_type = '{}', center_names = {}," \
                "original_p = {}, targeted_range = {})\n".format(
                    self.par_type, self.center_names,
                    self.original_p, self.targeted_range
                )
        elif self.par_type == "nbthole":
            return "DrudeParameter(par_type = '{}', center_names = {}," \
                "original_p = {}, targeted_range = {})\n".format(
                    self.par_type, self.center_names,
                    self.original_p, self.targeted_range
                )
        else:
            return "DrudeParameter with no change)\n"


def empty_drude_parameter():
    return DrudeParameter(
        lipidname='EMPTY', cgid=0, internal_id=0,
        par_type=None, center_names=[], original_p=0, targeted_range=[]
    )

# misc functions --------------------------------------------------------------


def get_sign(x):
    return math.copysign(1, x)


def gen_sensitivity_offset(paragroup, sign=1, percentage=1, lj_abs=False):
    """
    Function to generate the nonbonded parameter offset for a group of atoms.
    Args:
        paragroup: the parameter group generated by the
        fflip.chm.lipid.parse_nbgroup
        sign (int, optional): default to 1, the sign of the offset, can be 1/-1
        percentage (float): percentage% change of the offset,
        be careful and talk to Yalun before using
        lj_abs (bool): default to False, if use absolute value
        for LJ parameters (or percentage)
    Returns:
        the offset of the parameter (for the center atom(s))
    """
    pt = paragroup.par_type
    if pt == 'charge':
        num_atom = len(paragroup.center_names)
        return (0.01 * float(percentage) / num_atom) * sign * get_sign(
            paragroup.original_p
        )
    if pt == 'alpha':
        # return 0.0001 * float(percentage)  # 0.829 -> 0.929 if percentage is 1
        return 0.01 * float(percentage)
    if pt == 'thole' or pt == 'nbthole':
        # return 0.1 * float(percentage)
        return 0.01 * float(percentage)
    if pt == 'epsilon' or pt == 'sigma':
        if not lj_abs:
            return sign * 0.01 * float(percentage)
        else:
            if pt == 'epsilon':
                return sign * 0.005 * float(percentage) / paragroup.original_p
            else:
                return sign * 0.03 * float(percentage) / paragroup.original_p


def find_a_paragroup(groups, atom_name, par_type):
    succeed = False
    for g in groups:
        if g.par_type == par_type and atom_name in g.center_names:
            succeed = True
            group = g
    assert succeed
    return group

