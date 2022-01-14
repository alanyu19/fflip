# -*- coding: utf-8 -*-

import math


class NonbondedGroup(object):
    """
    A class that contains the info about the group for which one want to
    change the nobonded parameters but not the interaction form ...
    """
    def __init__(self, lipid_name, par_type, center_names, **kwargs):
        self.lipid = lipid_name
        self.par_type = par_type
        self.center_names = center_names
        if self.par_type == "charge":
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
            neib_names = []
            for neib in self.neighbors:
                neib_names.append(
                    neib + '_' + self.center_names[0] + '_' + self.par_type
                )
                name_list.append(
                    neib + '_' + self.center_names[0] + '_' + self.par_type
                )
            return [name, neib_names]
        else:
            return None

    def __repr__(self):
        # par_type, center_names,
        if self.par_type == "charge":
            return "NonbondedGroup(par_type='{}', center_names={}," \
                "original_p={}, targeted_range={}, neighbors={}," \
                " ron={})\n".format(
                    self.par_type, self.center_names, self.original_p,
                    self.targeted_range, self.neighbors, self.ron
                )
        elif self.par_type == "epsilon" or self.par_type == "sigma":
            return "NonbondedGroup(par_type='{}', center_names={}," \
                "original_p={}, targeted_range={})\n".format(
                    self.par_type, self.center_names,
                    self.original_p, self.targeted_range
                )
        else:
            return "nonbonded_group (no change)\n"


def empty_nonbonded_group():
    return NonbondedGroup(
        lipid_name='I am a fatty lipid', par_type=None, center_names=[], original_p=0, targeted_range=[]
    )


# -----------------------------------------------------------------------------
# DRUDE -----------------------------------------------------------------------


class DrudeParameter(object):
    def __init__(self, lipid_name, cgid, internal_id,
                 par_type, center_names, original_p, targeted_range,
                 neighbors=None, drude_particles=None, **kwargs):
        """
        Class contains
        Args:
            lipid_name: str, the lipid it belongs to
            cgid: int, the group id it belongs to (see dlipid)
            partype: can be 'charge, alpha, thole, sigma, epsilon
            center_names:
            neighbors:
            original_p:
            targeted_range:
            **kwargs:
        """
        self.lipid = lipid_name
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
            return "DrudeParameter(par_type='{}', center_names={}," \
                "original_p={}, targeted_range={})\n".format(
                    self.par_type, self.center_names,
                    self.original_p, self.targeted_range
                )
        elif self.par_type == "nbthole":
            return "DrudeParameter(par_type='{}', center_names={}," \
                "original_p={}, targeted_range={})\n".format(
                    self.par_type, self.center_names,
                    self.original_p, self.targeted_range
                )
        else:
            return "DrudeParameter with no change)\n"


def empty_drude_parameter():
    return DrudeParameter(
        lipid_name='EMPTY', cgid=0, internal_id=0,
        par_type=None, center_names=[], original_p=0, targeted_range=[]
    )

# misc functions --------------------------------------------------------------


