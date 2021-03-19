# -*- coding: utf-8 -*-


class LoadTrajectoryError(Exception):
    pass


class TrajNumberInsufficientError(LoadTrajectoryError):
    def __init__(self, ntraj, blksz):
        super(LoadTrajectoryError, self).__init__(
              "Number of short trajectories ({}) < block size ({}). Quitting.".
                  format(ntraj, blksz)
        )


class SetupSystemError(Exception):
    pass


class ParameterTypeNotExistentError(SetupSystemError):
    def __init__(self, wrong_par):
        super(SetupSystemError, self).__init__(
              "The parameter type given is '{}', "
              "but it did not match any of these 3: "
              "1. charge, 2. sigma, 3.epsilon ".format(wrong_par)
        )


class GlobalParameterForceTypeError(SetupSystemError):
    def __init__(self, wrong_par):
        super(SetupSystemError, self).__init__(
              "The force type given is '{}',"
              " but global parameter can only be set for nonbonded forces".
                  format(wrong_par)
        )


# ------------------------------ Reweighting ---------------------------------

class PropRewAvgPlotError(Exception):
    pass


class xyDimensionNotEqualError(PropRewAvgPlotError):
    def __init__(self, size1, size2):
        super(PropRewAvgPlotError, self).__init__(
              "X and Y axis have unequal dimension ({} vs. {}), Quitting.".
                  format(size1, size2)
        )


class GetPropOptionError(PropRewAvgPlotError):
    def __init__(self, option):
        super(PropRewAvgPlotError, self).__init__(
              "'{}' is not a good option. Supported option is "
              "1--plot only, 2--avg only, 3--both.".format(option)
        )


class StartEndTimeError(PropRewAvgPlotError):
    def __init__(self, start, end):
        super(PropRewAvgPlotError, self).__init__(
              "'{} ~ {} ns' is not a good range of simulation time.".format(
                  start, end
              )
        )


class ReweightingError(Exception):
    pass


class NoSimDataError(ReweightingError):
    def __init__(self):
        print('There is no simulated data, '
              'please run the reweight method first!')


class NoRewDataError(ReweightingError):
    def __init__(self):
        print('There is no  reweighted data, '
              'please run the reweight method first!')


class OtherReweightError(ReweightingError):
    def __int__(self, message):
        print(message)

# ------------------------ Other Coffe Module Related ------------------------


def omm_parse_fatal_error(filename):
    """Extract error message from OpenMM output.
    Arguments:
        file  --  a file containing OpenMM's stderr output
    """
    errmsg = ""
    in_msg = False
    with open(filename, "r") as f:
        for line in f:
            if line.strip().lower().startswith("fatal"):
                in_msg = True
            if in_msg:
                errmsg += line
    return errmsg


class OmmSimError(Exception):
    def __init__(self, exception, stderr_file):
        try:
            omm_msg = omm_parse_fatal_error(stderr_file)
            # assert omm_msg != ""
            msg = omm_msg + ("Details in {}\n\n\n".format(stderr_file))
            super(OmmSimError, self).__init__(msg)
        except Exception as e:
            super(OmmSimError, self).__init__(str(exception))
