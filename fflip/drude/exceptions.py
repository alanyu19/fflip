# -*- coding: utf-8 -*-


class GroupRepeatedError(Exception):
    def __init__(self, group):
        super().__init__("There was group {} already")
