#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys DC.
#
#  See LICENSE for licence details.

from typing import List, Optional

from hammer_vlsi import HammerSynthesisTool, HammerToolStep
from hammer_vlsi import SynopsysTool
from hammer_logging import HammerVLSILogging
import hammer_tech
from hammer_tech import HammerTechnologyUtils

import os
import re

class DC(HammerSynthesisTool, SynopsysTool):
    def fill_outputs(self) -> bool:
        # Check that the mapped.v exists if the synthesis run was successful
        # TODO: move this check upwards?
        mapped_v = "%s/results/%s.mapped.v" % (self.run_dir, self.top_module)
        if not os.path.isfile(mapped_v):
            raise ValueError("Output mapped verilog %s not found" % (mapped_v))  # better error?
        self.output_files = [mapped_v]
        return True

    @property
    def steps(self) -> List[HammerToolStep]:
        return self.make_steps_from_methods([
            self.main_step
        ])

    def tool_config_prefix(self) -> str:
        return "synthesis.dc"

    ## FIXME: this is abstract by SynopsysTool
    @property
    def post_synth_sdc(self) -> Optional[str]:
      return None

    # TODO(edwardw): move this to synopsys common
    def generate_tcl_preferred_routing_direction(self):
        """
        Generate a TCL fragment for setting preferred routing directions.
        """
        output = []

        # Suppress PSYN-882 ("Warning: Consecutive metal layers have the same preferred routing direction") while the layer routing is being built.
        output.append("set suppress_errors  [concat $suppress_errors  [list PSYN-882]]")

        # TODO This is broken with the API change to stackups ucb-bar/hammer#308
        #for library in self.technology.config.libraries:
        #    if library.metal_layers is not None:
        #        for layer in library.metal_layers:
        #            output.append("set_preferred_routing_direction -layers {{ {0} }} -direction {1}".format(layer.name,
                                                                                                            layer.preferred_routing_direction))

        output.append("set suppress_errors  [lminus $suppress_errors  [list PSYN-882]]")
        output.append("")  # Add newline at the end
        return "\n".join(output)

    def disable_congestion_map(self) -> None:
        """Disables the congestion map generation in rm_dc_scripts/dc.tcl since it requires a GUI and licences.
        """
        dc_tcl_path = os.path.join(self.run_dir, "rm_dc_scripts/dc.tcl")

        with open(dc_tcl_path) as f:
            dc_tcl = f.read()

        congestion_map_fragment = """
  # Use the following to generate and write out a congestion map from batch mode
  # This requires a GUI session to be temporarily opened and closed so a valid DISPLAY
  # must be set in your UNIX environment.

  if {[info exists env(DISPLAY)]} {
    gui_start
"""
        congestion_map_search = re.escape(congestion_map_fragment)
        # We want to capture & replace that condition.
        # Unfortunately, we can't replace within a group, so we'll have to replace around it.
        # e.g. foobarbaz -> foo123baz requires (foo)bar(baz) -> \1 123 \2 -> foo123baz
        cond = re.escape("[info exists env(DISPLAY)]")
        congestion_map_search_and_capture = "(" + congestion_map_search.replace(cond, ")(" + cond + ")(") + ")"

        output = re.sub(congestion_map_search_and_capture, "\g<1>false\g<3>", dc_tcl)

        f = open(dc_tcl_path, 'w')
        f.write(output)
        f.close()

    def main_step(self) -> bool:
        # TODO(edwardw): move most of this to Synopsys common since it's not DC-specific.
        # Locate reference methodology tarball.
        synopsys_rm_tarball = self.get_synopsys_rm_tarball("DC")

        # Locate DC binary.
        dc_bin = self.get_setting("synthesis.dc.dc_bin")
        if not os.path.isfile(dc_bin):
            self.logger.error("DC binary not found as expected at {0}".format(dc_bin))
            # TODO(edwardw): think about how to pass extra/more useful info like the specific type of error (e.g. FileNotFoundError)
            return False

        # Load input files.
        if not self.check_input_files([".v", ".sv"]):
            return False

        input_files = list(self.input_files)  # type: List[str]
        # Add any verilog_synth wrappers (which are needed in some technologies e.g. for SRAMs) which need to be
        # synthesized.
        input_files += self.technology.read_libs([
            hammer_tech.filters.verilog_synth_filter
        ], HammerTechnologyUtils.to_plain_item)

        # Generate preferred_routing_directions.
        preferred_routing_directions_fragment = os.path.join(self.run_dir, "preferred_routing_directions.tcl")
        with open(preferred_routing_directions_fragment, "w") as f:
            f.write(self.generate_tcl_preferred_routing_direction())

        # Generate clock constraints.
        clock_constraints_fragment = os.path.join(self.run_dir, "clock_constraints_fragment.tcl")
        with open(clock_constraints_fragment, "w") as f:
            f.write(self.sdc_clock_constraints)

        # Get libraries.
        lib_args = self.technology.read_libs([
            hammer_tech.filters.timing_db_filter._replace(tag="lib"),
            hammer_tech.filters.milkyway_lib_dir_filter._replace(tag="milkyway"),
            hammer_tech.filters.tlu_max_cap_filter._replace(tag="tlu_max"),
            hammer_tech.filters.tlu_min_cap_filter._replace(tag="tlu_min"),
            hammer_tech.filters.tlu_map_file_filter._replace(tag="tlu_map"),
            hammer_tech.filters.milkyway_techfile_filter._replace(tag="tf")
        ], HammerTechnologyUtils.to_command_line_args)

        # Pre-extract the tarball (so that we can make TCL modifications in Python)
        self.run_executable([
            "tar", "-xf", synopsys_rm_tarball, "-C", self.run_dir, "--strip-components=1"
        ])

        # Disable the DC congestion map if needed.
        if not self.get_setting("synthesis.dc.enable_congestion_map"):
            self.disable_congestion_map()

        compile_args = self.get_setting("synthesis.dc.compile_args")
        if not compile_args:
          compile_args = []

        # Build args.
        args = [
            os.path.join(self.tool_dir, "tools", "run-synthesis"),
            "--dc", dc_bin,
            "--clock_constraints_fragment", clock_constraints_fragment,
            "--preferred_routing_directions_fragment", preferred_routing_directions_fragment,
            "--find_regs_tcl", os.path.join(self.tool_dir, "tools", "find-regs.tcl"),
            "--run_dir", self.run_dir,
            "--top", self.top_module
        ]
        args.extend(input_files)  # We asserted these are Verilog above
        args.extend(lib_args)
        for compile_arg in compile_args:
          args.extend(["--compile_arg", compile_arg])

        # Temporarily disable colours/tag to make DC run output more readable.
        # TODO: think of a more elegant way to do this?
        HammerVLSILogging.enable_colour = False
        HammerVLSILogging.enable_tag = False
        self.run_executable(args)  # TODO: check for errors and deal with them
        HammerVLSILogging.enable_colour = True
        HammerVLSILogging.enable_tag = True
        return True


tool = DC
