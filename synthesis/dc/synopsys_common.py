#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys Tools.
#
#  See LICENSE for licence details.

from typing import List

import os

import hammer_tech
from hammer_tech import HammerTechnologyUtils
from hammer_vlsi import SynopsysTool

#TODO: this probably needs to be just SynopsysTool
class SynopsysCommon(SynopsysTool):
    @property
    def script_dir(self) -> str:
        dirname = os.path.join(self.run_dir, "scripts")
        os.makedirs(dirname, exist_ok=True)
        return dirname

    @property
    def report_dir(self) -> str:
        dirname = os.path.join(self.run_dir, "reports")
        os.makedirs(dirname, exist_ok=True)
        return dirname

    @property
    def result_dir(self) -> str:
        dirname = os.path.join(self.run_dir, "results")
        os.makedirs(dirname, exist_ok=True)
        return dirname

    @property
    def timing_dbs(self) -> List[str]:
        # Gather/load libraries.
        return self.technology.read_libs(
            [hammer_tech.filters.timing_db_filter],
            HammerTechnologyUtils.to_plain_item)

    @property
    def milkyway_lib_dirs(self) -> List[str]:
        return self.technology.read_libs(
            [hammer_tech.filters.milkyway_lib_dir_filter],
            HammerTechnologyUtils.to_plain_item)

    @property
    def milkyway_techfiles(self) -> List[str]:
        return self.technology.read_libs(
            [hammer_tech.filters.milkyway_techfile_filter],
            HammerTechnologyUtils.to_plain_item)

    @property
    def tlu_max_caps(self) -> List[str]:
        return self.technology.read_libs(
            [hammer_tech.filters.tlu_max_cap_filter],
            HammerTechnologyUtils.to_plain_item)

    @property
    def tlu_min_caps(self) -> List[str]:
        return self.technology.read_libs(
            [hammer_tech.filters.tlu_min_cap_filter],
            HammerTechnologyUtils.to_plain_item)

    @property
    def tlu_map(self) -> List[str]:
        return self.technology.read_libs(
            [hammer_tech.filters.tlu_map_file_filter],
            HammerTechnologyUtils.to_plain_item)

    @property
    def verilog(self) -> List[str]:
        return [v for v in list(self.input_files) if v.endswith(".v")]
