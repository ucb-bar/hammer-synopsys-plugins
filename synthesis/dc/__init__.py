#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys DC.
#
#  See LICENSE for licence details.

from typing import List, Optional, Dict

import os
import re

from hammer_vlsi import HammerSynthesisTool, HammerToolStep
from hammer_logging import HammerVLSILogging
import hammer_tech
from hammer_tech import HammerTechnologyUtils
from .synopsys_common import SynopsysCommon

class DC(HammerSynthesisTool, SynopsysCommon):
    def fill_outputs(self) -> bool:
        # Check that the mapped.v exists if the synthesis run was successful
        # TODO: move this check upwards?
        mapped_v = os.path.join(self.result_dir, self.top_module + ".mapped.v")
        if not os.path.isfile(mapped_v):
            raise ValueError("Output mapped verilog %s not found" % (mapped_v))  # better error?
        self.output_files = [mapped_v]
        return True

    def tool_config_prefix(self) -> str:
        return "synthesis.dc"

    @property
    def post_synth_sdc(self) -> Optional[str]:
        return os.path.join(self.result_dir, self.top_module + ".mapped.sdc")

    @property
    def steps(self) -> List[HammerToolStep]:
        return self.make_steps_from_methods([
            self.init_environment,
            self.elaborate_design,
            self.apply_constraints,
            self.optimize_design,
            self.generate_reports,
            self.write_outputs,
        ])

    def do_post_steps(self) -> bool:
        assert super().do_post_steps()
        return self.run_design_compiler()

    @property
    def output(self) -> List[str]:
        """
        Buffered output to be put into dc.tcl.
        """
        return self.attr_getter("_output", [])

    def append(self, cmd: str) -> None:
        self.tcl_append(cmd, self.output)


    def init_environment(self) -> bool:
        # The following setting removes new variable info messages from the end of the log file
        self.append("set_app_var sh_new_variable_message false")

        # Actually use specified number of cores
        self.append("set disable_multicore_resource_checks true")
        self.append("set_host_options -max_cores %d" % self.get_setting("vlsi.core.max_threads"))

        # Change alib_library_analysis_path to point to a central cache of analyzed libraries
        # to save runtime and disk space.  The following setting only reflects the
        # default value and should be changed to a central location for best results.
        self.append("set_app_var alib_library_analysis_path alib")

        # Search Path Setup
        self.append("set_app_var search_path \". %s $search_path\"" % self.result_dir)

        # Library setup
        for db in self.timing_dbs:
            if not os.path.exists(db):
                self.logger.error("Cannot find %s" % db)
                return False
        self.append("set_app_var target_library \"%s\"" % ' '.join(self.timing_dbs))
        self.append("set_app_var synthetic_library dw_foundation.sldb")
        self.append("set_app_var link_library \"* $target_library $synthetic_library\"")

        # For designs that don't have tight QoR constraints and don't have register retiming,
        # you can use the following variable to enable the highest productivity single pass flow.
        # This flow modifies the optimizations to make verification easier.
        # This variable setting should be applied prior to reading in the RTL for the design.
        self.append("set_app_var simplified_verification_mode false")
        self.append("set_svf results/%s.mapped.svf" % self.top_module)

        return True

    def elaborate_design(self) -> bool:
        # Add any verilog_synth wrappers
        # (which are needed in some technologies e.g. for SRAMs)
        # which need to be synthesized.
        verilog = self.verilog + self.technology.read_libs([
            hammer_tech.filters.verilog_synth_filter
        ], HammerTechnologyUtils.to_plain_item)
        for v in verilog:
            if not os.path.exists(v):
                self.logger.error("Cannot find %s" % v)
                return False
        # Read RTL
        self.append("define_design_lib WORK -path ./WORK")
        self.append("analyze -format verilog \"%s\"" % ' '.join(verilog))

        # Elabrate design
        self.append("elaborate %s" % self.top_module)
        return True

    def apply_constraints(self) -> bool:
        # Generate clock
        clocks = [clock.name for clock in self.get_clock_ports()]
        self.append(self.sdc_clock_constraints)

        # Set ungroup
        for module in self.get_setting('vlsi.inputs.no_ungroup'):
            self.append("set_ungroup [get_designs %s] false" % module)

        # Set retmining
        for module in self.get_setting("vlsi.inputs.retimed_modules"):
            self.append(' '.join([
                "set_optimize_registers", "true",
                "-design", module,
                "-clock", "{%s}" % ' '.join(clocks)
            ] + self.get_setting("synthesis.dc.retiming_args")))

        # Create Default Path Groups
        self.append("""
set ports_clock_root [filter_collection [get_attribute [get_clocks] sources] object_class==port]
group_path -name REGOUT -to [all_outputs]
group_path -name REGIN -from [remove_from_collection [all_inputs] ${ports_clock_root}]
group_path -name FEEDTHROUGH -from [remove_from_collection [all_inputs] ${ports_clock_root}] -to [all_outputs]
""")
        # Prevent assignment statements in the Verilog netlist.
        self.append("set_fix_multiple_port_nets -all -buffer_constants")

        return True


    def optimize_design(self) -> bool:
        # Optimize design
        self.append("compile_ultra %s" % ' '.join(self.get_setting("synthesis.dc.compile_args")))
        self.append("change_names -rules verilog -hierarchy")
        # Write and close SVF file and make it available for immediate use
        self.append("set_svf -off")
        return True

    def generate_reports(self) -> bool:
        self.append("""
report_reference -hierarchy > \\
    {report_dir}/{design_name}.mapped.report_reference.out
report_qor > \\
    {report_dir}/{design_name}.mapped.qor.rpt
report_area -nosplit > \\
    {report_dir}/{design_name}.mapped.area.rpt
report_timing -max_paths 500 -nworst 10 -input_pins -capacitance \\
    -significant_digits 4 -transition_time -nets -attributes -nosplit > \\
    {report_dir}/{design_name}.mapped.timing.rpt

report_power -nosplit > \\
    {report_dir}/{design_name}.mapped.power.rpt
report_clock_gating -nosplit > \\
    {report_dir}/{design_name}.mapped.clock_gating.rpt
""".format(report_dir=self.report_dir, design_name=self.top_module))
        return True

    def write_outputs(self) -> bool:
        self.append("""
write -format verilog -hierarchy -output \\
    {result_dir}/{design_name}.mapped.v
write -format ddc -hierarchy -output \\
    {result_dir}/{design_name}.mapped.ddc
write_sdc -nosplit \\
    {result_dir}/{design_name}.mapped.sdc
""".format(result_dir=self.result_dir, design_name=self.top_module))
        return True

    @property
    def env_vars(self) -> Dict[str, str]:
        env = dict(super().env_vars)
        env["PATH"] = "%s:%s" % (
            os.path.dirname(self.get_setting("synthesis.dc.dc_bin")),
            os.environ["PATH"])
        return env

    def run_design_compiler(self) -> bool:
        HammerVLSILogging.enable_colour = False
        HammerVLSILogging.enable_tag = False
        dc_bin = os.path.basename(self.get_setting("synthesis.dc.dc_bin"))
        dc_tcl = os.path.join(self.script_dir, "dc.tcl")
        with open(dc_tcl, 'w') as _f:
            _f.write('\n'.join(self.output))
            _f.write('\nexit')
        args = [dc_bin, "-64bit", "-f", dc_tcl]
        # TODO: check outputs from lines?
        lines = self.run_executable(args, self.run_dir)
        HammerVLSILogging.enable_colour = True
        HammerVLSILogging.enable_tag = True
        return True

tool = DC
