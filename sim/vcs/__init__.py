#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys VCS
#
#  See LICENSE for license details.

from hammer_vlsi import HammerSimTool, HammerToolStep
from hammer_vlsi import SynopsysTool
from hammer_logging import HammerVLSILogging

from typing import Dict, List, Optional, Callable, Tuple

from hammer_vlsi import SimulationLevel

import hammer_utils
import hammer_tech
from hammer_tech import HammerTechnologyUtils

import os
import re
import shutil

class VCS(HammerSimTool, SynopsysTool):

    def post_synth_sdc(self) -> Optional[str]:
        pass

    def tool_config_prefix(self) -> str:
        return "sim.vcs"

    @property
    def steps(self) -> List[HammerToolStep]:
        return self.make_steps_from_methods([
            self.write_gl_files,
            self.run_vcs,
            self.run_simulation
            ])

    def benchmark_run_dir(self, bmark_path: str) -> str:
        """Generate a benchmark run directory."""
        # TODO(ucb-bar/hammer#462) this method should be passed the name of the bmark rather than its path
        bmark = os.path.basename(bmark_path)
        return os.path.join(self.run_dir, bmark)

    @property
    def force_regs_file_path(self) -> str:
        return os.path.join(self.run_dir, "force_regs.ucli")

    @property
    def access_tab_file_path(self) -> str:
        return os.path.join(self.run_dir, "access.tab")

    @property
    def simulator_executable_path(self) -> str:
        return os.path.join(self.run_dir, "simv")

    @property
    def run_tcl_path(self) -> str:
        return os.path.join(self.run_dir, "run.tcl")

    @property
    def env_vars(self) -> Dict[str, str]:
        v = dict(super().env_vars)
        v["VCS_HOME"] = self.get_setting("sim.vcs.vcs_home")
        v["SNPSLMD_LICENSE_FILE"] = self.get_setting("synopsys.SNPSLMD_LICENSE_FILE")
        return v

    def get_verilog_models(self) -> List[str]:
        verilog_sim_files = self.technology.read_libs([
            hammer_tech.filters.verilog_sim_filter
        ], hammer_tech.HammerTechnologyUtils.to_plain_item)
        return verilog_sim_files

    def write_gl_files(self) -> bool:
        if self.level == SimulationLevel.RTL:
            return True

        tb_prefix = self.get_setting("sim.inputs.tb_dut")
        force_val = self.get_setting("sim.inputs.gl_register_force_value")

        seq_cells = self.seq_cells

        with open(self.access_tab_file_path, "w") as f:
            for cell in seq_cells:
                f.write("acc=wn:{cell_name}\n".format(cell_name=cell))

        all_regs = self.all_regs

        with open(self.force_regs_file_path, "w") as f:
            for reg in all_regs:
                path = reg["path"]
                path = '.'.join(path.split('/'))
                pin = reg["pin"]
                f.write("force -deposit {" + tb_prefix + "." + path + " ." + pin + "} " + str(force_val) + "\n")

        return True

    def run_vcs(self) -> bool:
        # run through inputs and append to CL arguments
        vcs_bin = self.get_setting("sim.vcs.vcs_bin")
        if not os.path.isfile(vcs_bin):
          self.logger.error("VCS binary not found as expected at {0}".format(vcs_bin))
          return False

        if not self.check_input_files([".v", ".sv", ".so", ".cc", ".c"]):
          return False

        top_module = self.top_module
        compiler_opts = self.get_setting("sim.inputs.compiler_opts", [])
        # TODO(johnwright) sanity check the timescale string
        timescale = self.get_setting("sim.inputs.timescale")
        input_files = list(self.input_files)
        options = self.get_setting("sim.inputs.options", [])
        defines = self.get_setting("sim.inputs.defines", [])
        access_tab_filename = self.access_tab_file_path
        tb_name = self.get_setting("sim.inputs.tb_name")

        # Build args
        args = [
          vcs_bin,
          "-full64"
        ]

        if timescale is not None:
            args.append('-timescale={}'.format(timescale))

        # Add in options we pass to the C++ compiler
        args.extend(['-CC', '-I$(VCS_HOME)/include'])
        for compiler_opt in compiler_opts:
            args.extend(['-CC', compiler_opt])

        # black box options
        args.extend(options)

        # Add in all input files
        args.extend(input_files)

        # Note: we always want to get the verilog models because most real designs will instantate a few
        # tech-specific cells in the source RTL (IO cells, clock gaters, etc.)
        args.extend(self.get_verilog_models())

        for define in defines:
            args.extend(['+define+' + define])

        if self.level == SimulationLevel.GateLevel:
            args.extend(['-P'])
            args.extend([access_tab_filename])
            args.extend(['-debug'])
            if self.get_setting("sim.inputs.timing_annotated"):
                args.extend(["+neg_tchk"])
                args.extend(["+sdfverbose"])
                args.extend(["-negdelay"])
                args.extend(["-sdf"])
                args.extend(["max:{top}:{sdf}".format(run_dir=self.run_dir, top=top_module, sdf=self.sdf_file)])
            else:
                args.extend(["+notimingcheck"])
                args.extend(["+delay_mode_zero"])

        args.extend(["-top", tb_name])

        args.extend(['-o', self.simulator_executable_path])

        HammerVLSILogging.enable_colour = False
        HammerVLSILogging.enable_tag = False

        # Delete an old copy of the simulator if it exists
        if os.path.exists(self.simulator_executable_path):
            os.remove(self.simulator_executable_path)

        # Remove the csrc directory (otherwise the simulator will be stale)
        if os.path.exists(os.path.join(self.run_dir, "csrc")):
            shutil.rmtree(os.path.join(self.run_dir, "csrc"))

        # Generate a simulator
        self.run_executable(args, cwd=self.run_dir)

        HammerVLSILogging.enable_colour = True
        HammerVLSILogging.enable_tag = True

        return os.path.exists(self.simulator_executable_path)

    def run_simulation(self) -> bool:
        if not self.get_setting("sim.inputs.execute_sim"):
            self.logger.warning("Not running any simulations because sim.inputs.execute_sim is unset.")
            return True

        top_module = self.top_module
        exec_flags_prepend = self.get_setting("sim.inputs.execution_flags_prepend", [])
        exec_flags = self.get_setting("sim.inputs.execution_flags", [])
        exec_flags_append = self.get_setting("sim.inputs.execution_flags_append", [])
        force_regs_filename = self.force_regs_file_path

        if self.level == SimulationLevel.GateLevel:
            with open(self.run_tcl_path, "w") as f:
                find_regs_run_tcl = []
                find_regs_run_tcl.append("source " + force_regs_filename)
                find_regs_run_tcl.append("run")
                find_regs_run_tcl.append("exit")
                f.write("\n".join(find_regs_run_tcl))

        for benchmark in self.benchmarks:
            if not os.path.isfile(benchmark):
              self.logger.error("benchmark not found as expected at {0}".format(vcs_bin))
              return False

        # setup simulation arguments
        args = [ self.simulator_executable_path ]
        args.extend(exec_flags_prepend)
        args.extend(exec_flags)
        if self.level == SimulationLevel.GateLevel:
            args.extend(["-ucli", "-do", self.run_tcl_path])
        args.extend(exec_flags_append)

        HammerVLSILogging.enable_colour = False
        HammerVLSILogging.enable_tag = False

        # TODO(johnwright) We should optionally parallelize this in the future.
        for benchmark in self.benchmarks:
            bmark_run_dir = self.benchmark_run_dir(benchmark)
            # Make the rundir if it does not exist
            hammer_utils.mkdir_p(bmark_run_dir)
            self.run_executable(args + [benchmark], cwd=bmark_run_dir)

        if self.benchmarks == []:
            self.run_executable(args, cwd=self.run_dir)

        HammerVLSILogging.enable_colour = True
        HammerVLSILogging.enable_tag = True

        return True

tool = VCS
