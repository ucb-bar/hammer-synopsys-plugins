#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys DC.
#
#  See LICENSE for licence details.

from typing import List, Optional, Dict, Any

import os
import re

from hammer.vlsi import HammerSynthesisTool, HammerToolStep
from hammer.logging import HammerVLSILogging
import hammer.tech
from hammer.tech import HammerTechnologyUtils
from .synopsys_common import SynopsysCommon

class DC(HammerSynthesisTool, SynopsysCommon):
    def fill_outputs(self) -> bool:
        # Check that the regs paths were written properly if the write_regs step was run
        self.output_seq_cells = self.all_cells_path
        self.output_all_regs = self.all_regs_path
        if self.ran_write_regs:
            if not os.path.isfile(self.all_cells_path):
                raise ValueError("Output find_regs_cells.json %s not found" % (self.all_cells_path))

            if not os.path.isfile(self.all_regs_path):
                raise ValueError("Output find_regs_paths.json %s not found" % (self.all_regs_path))

            if not self.process_reg_paths(self.all_regs_path):
                self.logger.error("Failed to process all register paths")
        else:
            self.logger.info("Did not run write_regs")
        # Check that the mapped.v exists if the synthesis run was successful
        # TODO: move this check upwards?
        mapped_v = os.path.join(self.result_dir, self.top_module + ".mapped.v")
        if not os.path.isfile(mapped_v):
            raise ValueError("Output mapped verilog %s not found" % (mapped_v))  # better error?
        self.output_files = [mapped_v]
        # DC does not 
        self.output_sdc = self.post_synth_sdc
        self.sdf_file = self.output_sdf_path
        if self.ran_write_outputs:
            if not os.path.isfile(mapped_v):
                raise ValueError("Output mapped verilog %s not found" % (mapped_v)) # better error?

            if not os.path.isfile(self.output_sdc):
                self.logger.warning("Output SDC %s not found" % (self.mapped_sdc_path)) # better error?

            if not os.path.isfile(self.output_sdf_path):
                self.logger.warning("Output SDF %s not found" % (self.output_sdf_path))
        else:
            self.logger.info("Did not run write_outputs")

        return True

    def tool_config_prefix(self) -> str:
        return "synthesis.dc"

    @property
    def all_regs_path(self) -> str:
        return os.path.join(self.run_dir, "find_regs_paths.json")

    @property
    def all_cells_path(self) -> str:
        return os.path.join(self.run_dir, "find_regs_cells.json")

    @property
    def ran_write_regs(self) -> bool:
        """The write_regs step sets this to True if it was run."""
        return self.attr_getter("_ran_write_regs", False)

    @ran_write_regs.setter
    def ran_write_regs(self, val: bool) -> None:
        self.attr_setter("_ran_write_regs", val)
    
    @property
    def ran_write_outputs(self) -> bool:
        """The write_ouputs step sets this to True if it was run."""
        return self.attr_getter("_ran_write_outputs", False)

    @ran_write_outputs.setter
    def ran_write_outputs(self, val: bool) -> None:
        self.attr_setter("_ran_write_outputs", val)

    def export_config_outputs(self) -> Dict[str, Any]:
        outputs = dict(super().export_config_outputs())
        # TODO(edwardw): find a "safer" way of passing around these settings keys.
        outputs["synthesis.outputs.sdc"] = self.output_sdc
        outputs["synthesis.outputs.seq_cells"] = self.output_seq_cells
        outputs["synthesis.outputs.all_regs"] = self.output_all_regs
        outputs["synthesis.outputs.sdf_file"] = self.output_sdf_path
        return outputs
    
    @property
    def post_synth_sdc(self) -> Optional[str]:
        return os.path.join(self.result_dir, self.top_module + ".mapped.sdc")
    
    @property
    def output_sdf_path(self) -> str:
        return os.path.join(self.run_dir, "{top}.mapped.sdf".format(top=self.top_module)) 

    @property
    def steps(self) -> List[HammerToolStep]:
        return self.make_steps_from_methods([
            self.init_environment,
            self.elaborate_design,
            self.apply_constraints,
            self.insert_dft,
            self.optimize_design,
            self.generate_reports,
            self.generate_dft_reports,
            self.write_outputs,
            self.write_regs,
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
            hammer.tech.filters.verilog_synth_filter
        ], HammerTechnologyUtils.to_plain_item)
        for v in verilog:
            if not os.path.exists(v):
                self.logger.error("Cannot find %s" % v)
                return False
        # Read RTL
        self.append("define_design_lib WORK -path ./WORK")
        self.append("analyze -format sverilog \"%s\"" % ' '.join(verilog))

        # Elaborate design
        self.append("elaborate %s" % self.top_module)
        
        # Se the current design
        self.append("current_design %s" % self.top_module)

        # Link the design for possible unresolved errors
        self.append("link")

        # Set Rams as black boxes
        #self.append("foreach_in_collection ram [get_designs *ram*] \{ ")
        #self.append("set_attribute $ram is_black_box true")
        #self.append("\} ")
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
        self.ran_write_outputs = True 
        return True

    def generate_dft_reports(self) -> bool:
        self.append("""
write_test_protocol -output {result_dir}/{design_name}_test_protocol.spf
""".format(result_dir=self.result_dir, design_name=self.top_module))
        self.append("""
write_scan_def -output {result_dir}/{design_name}_report_dft.scandef
""".format(result_dir=self.result_dir, design_name=self.top_module))
        
        return True

    def write_regs(self) -> bool:
        """write regs info to be read in for simulation register forcing"""
        if self.hierarchical_mode.is_nonleaf_hierarchical():
            self.append(self.child_modules_tcl())
        self.append(self.write_regs_tcl())
        self.ran_write_regs = True
        return True

    def insert_dft(self) -> bool:
        clocks = [clock.name for clock in self.get_clock_ports()]
        resets = [reset.name for reset in self.get_reset_ports()]
        reset_active_negated = [0 if reset.active_negated else 1 for reset in self.get_reset_ports()]
        self.append("set compile_timing_high_effort true")
        self.append("set compile_delete_unloaded_sequential_cells true")
        self.append("set_scan_configuration -style multiplexed_flip_flop")
        self.append("compile -scan  -gate_clock -area_effor medium -map_effort medium")
        self.append("set_scan_configuration -chain_count 1  -create_test_clocks_by_system_clock_domain true")
        self.append("""
set_dft_signal -view existing_dft -type ScanClock -port {clock}  -timing [list 45 95] -active_state 1 -connect_to {clock}
""".format(clock=clocks[0]))
        self.append("create_port test_si -direction in")
        self.append("create_port test_se -direction in")
        self.append("create_port test_so -direction out")
        self.append("""
set_dft_signal -view spec -type Reset -port {reset} -active_state {state}
""".format(reset=resets[0],state=reset_active_negated[0]))
        self.append("set_dft_signal -view spec -type ScanDataIn -port test_si ")
        self.append("set_dft_signal -view spec -type ScanDataOut -port test_so")
        self.append("set_dft_signal -view spec -type ScanEnable -port test_se -active_state 1")
        self.append("create_test_protocol")
        self.append("dft_drc -verbose")
        self.append("preview_dft")
        self.append("insert_dft")
        self.append("check_scan")
        self.append("dft_drc")
        self.append("check_design")
        
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