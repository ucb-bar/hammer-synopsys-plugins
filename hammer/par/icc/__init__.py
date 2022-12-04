#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys ICC.
#
#  See LICENSE for licence details.

import os
import sys
import errno
from shutil import copyfile
from os.path import dirname
from typing import List, Optional,Callable, Tuple, Set, Any, cast, Dict
from decimal import Decimal

from hammer.utils import get_or_else, optional_map
import hammer.tech as hammer_tech
from hammer.tech import RoutingDirection, Metal, LibraryFilter
from hammer.vlsi import HammerTool, HammerToolHookAction, HammerPlaceAndRouteTool, HammerToolStep, PlacementConstraintType, HierarchicalMode, ObstructionType, Margins, Supply, PlacementConstraint, MMMCCornerType
from hammer.vlsi import SynopsysTool
from hammer.logging import HammerVLSILogging
from hammer.tech import HammerTechnologyUtils
from hammer.tech.specialcells import CellType, SpecialCell
from hammer.vlsi.constraints import BumpsDefinition

class ICC(HammerPlaceAndRouteTool, SynopsysTool):

    def export_config_outputs(self) -> Dict[str, Any]:
        outputs = dict(super().export_config_outputs())
        outputs["par.outputs.gds_file"] = self.output_gds_filename
        outputs["par.outputs.netlist"] = self.output_netlist_filename
        outputs["par.outputs.sim_netlist"] = self.output_sim_netlist_filename
        outputs["par.outputs.def_file"] = self.output_def_path
        outputs["par.outputs.sdf_file"] = self.output_sdf_path
        outputs["par.outputs.spefs"] = self.output_spef_paths
        return outputs

    def fill_outputs(self) -> bool:

        # TODO: To implement ILM run - needs to be included after implementation of hierarchical flow

        # Check for list of par outputs
        self.output_gds = self.output_gds_filename
        self.output_netlist = self.output_netlist_filename
        self.output_sim_netlist = self.output_sim_netlist_filename
        self.sdf_file = self.output_sdf_path
        self.spef_files = self.output_spef_paths

        if self.ran_write_design:
            if not os.path.isfile(self.output_gds_filename):
                raise ValueError("Output GDS %s not found" % (self.output_gds_filename))

            if not os.path.isfile(self.output_netlist_filename):
                raise ValueError("Output netlist %s not found" % (self.output_netlist_filename))

            if not os.path.isfile(self.output_sim_netlist_filename):
                raise ValueError("Output sim netlist %s not found" % (self.output_sim_netlist_filename))

            if not os.path.isfile(self.output_sdf_path):
                raise ValueError("Output SDF %s not found" % (self.output_sdf_path))

            for spef_path in self.output_spef_paths:
                if not os.path.isfile(spef_path):
                    raise ValueError("Output SPEF %s not found" % (spef_path))
        else:
            self.logger.info("Did not run write_design")

        return True

    @property
    def output_gds_filename(self) -> str:
        return os.path.join(self.run_dir, "{top}.gds".format(top=self.top_module))

    @property
    def output_netlist_filename(self) -> str:
        return os.path.join(self.run_dir, "{top}.lvs.v".format(top=self.top_module))

    @property
    def output_sim_netlist_filename(self) -> str:
        return os.path.join(self.run_dir, "{top}.sim.v".format(top=self.top_module)) 

    @property
    def output_sdf_path(self) -> str:
        return os.path.join(self.run_dir, "{top}.par.sdf".format(top=self.top_module))

    @property
    def output_def_path(self) -> str:
        return os.path.join(self.run_dir, "{top}.par.def".format(top=self.top_module))   

    @property
    def output_spef_paths(self) -> List[str]:
            corners = self.get_mmmc_corners()
            # TODO: To add corner specific paths after implementing MCMM flow
            spef_paths = []
            spef_paths.append(os.path.join(self.run_dir, "{top}.par.spef.max".format(top=self.top_module)))
            spef_paths.append(os.path.join(self.run_dir, "{top}.par.spef.min".format(top=self.top_module)))
            return spef_paths

    @property
    def env_vars(self) -> Dict[str, str]:
        v = dict(super().env_vars)
        v["ICC_BIN"] = self.get_setting("par.icc.icc_bin")
        return v    

    @property
    def _step_transitions(self) -> List[Tuple[str, str]]:
        """
        Private helper property to keep track of which steps we ran so that we
        can create symlinks.
        This is a list of (pre, post) steps
        """
        return self.attr_getter("__step_transitions", [])

    @_step_transitions.setter
    def _step_transitions(self, value: List[Tuple[str, str]]) -> None:
        self.attr_setter("__step_transitions", value)

    def do_pre_steps(self, first_step: HammerToolStep) -> bool:
        assert super().do_pre_steps(first_step)
        # Restore from the last checkpoint if we're not starting over.
        if first_step.name != "init_design":
           self.verbose_append("open_mw_lib {top}_LIB".format(top=self.top_module)) 
           self.verbose_append("open_mw_cel pre_{step}".format(step=first_step.name)) 
        return True

    def do_between_steps(self, prev: HammerToolStep, next: HammerToolStep) -> bool:
        assert super().do_between_steps(prev, next)
        if (prev.name != "init_design"):
            # Execute the pg_connection commands
            self.verbose_append("\n".join(self.pg_connection())) 
        # Write a checkpoint to disk.
        self.verbose_append("save_mw_cel -as pre_{step}".format(step=next.name))
        self._step_transitions = self._step_transitions + [(prev.name, next.name)]
        return True

    def do_post_steps(self) -> bool:
        assert super().do_post_steps()
        # TODO: this doesn't work if you're only running the very last step
        if len(self._step_transitions) > 0:
            last = "post_{step}".format(step=self._step_transitions[-1][1])
            self.verbose_append("save_mw_cel -as {last}".format(last=last))

        return self.run_icc()

    def get_tool_hooks(self) -> List[HammerToolHookAction]:
        return [self.make_persistent_hook(icc_global_settings)]

    @property
    def steps(self) -> List[HammerToolStep]:
        steps = [
            self.init_design,
            self.place_pins,
            self.floorplan_design,
            self.place_bumps,
            self.place_tap_cells,
            self.power_straps,
            self.place_opt_design,
            self.clock_tree,
            self.cts_opt,
            self.route_design,
            self.add_fillers
        ]
        write_design_step = [
            self.write_design
        ]  # type: List[Callable[[], bool]]             
        
        # Future TODO: Add hierarchical mode
        if self.hierarchical_mode == HierarchicalMode.Flat:
            # Nothing to do
            pass
        else:
            raise NotImplementedError("HierarchicalMode not implemented: " + str(self.hierarchical_mode))
        
        return self.make_steps_from_methods(steps + write_design_step)

    def tool_config_prefix(self) -> str:
        return "par.icc"

    def init_design(self) -> bool:
        """Initialize the design."""
        verbose_append = self.verbose_append

        # Common settings for the run
        # Additional common flags can be added if needed
        verbose_append("set hdlin_enable_rtldrc_info true")
        verbose_append("set enable_recovery_removal_arcs true")

        # Gather/load libraries.
        timing_dbs = ' '.join(self.technology.read_libs([hammer_tech.filters.timing_db_filter], HammerTechnologyUtils.to_plain_item))
        milkyway_lib_dirs = ' '.join(self.technology.read_libs([hammer_tech.filters.milkyway_lib_dir_filter],
HammerTechnologyUtils.to_plain_item))
        milkyway_techfiles = ' '.join(self.technology.read_libs([hammer_tech.filters.milkyway_techfile_filter], HammerTechnologyUtils.to_plain_item))
        tlu_max_caps = ' '.join(self.technology.read_libs([hammer_tech.filters.tlu_max_cap_filter], HammerTechnologyUtils.to_plain_item))
        tlu_min_caps = ' '.join(self.technology.read_libs([hammer_tech.filters.tlu_min_cap_filter], HammerTechnologyUtils.to_plain_item))
        tlu_map = ' '.join(self.technology.read_libs([hammer_tech.filters.tlu_map_file_filter], HammerTechnologyUtils.to_plain_item))       
 
        if timing_dbs == "":
            self.logger.error("No timing dbs (libs) specified!")
            return False
        if milkyway_lib_dirs == "":
            self.logger.error("No milkyway lib dirs specified!")
            return False
        if milkyway_techfiles == "":
            self.logger.error("No milkyway tech files specified!")
            return False
        if tlu_max_caps == "" and tlu_min_caps == "":
            self.logger.error("No tlu+ cap files specified!")
            return False

        # Load input files.
        if not self.check_input_files([".v"]):
            return False
        
        verilog_args = list(self.input_files)
        if len(verilog_args) < 1:
            self.logger.error("No post-synthesis Verilog netlists specified!")
            return False

        # Create Milkyway Library
        verbose_append("create_mw_lib -technology {techfile} -mw_reference_library {mw_ref} {top}_LIB".format(
                        techfile = milkyway_techfiles, mw_ref = milkyway_lib_dirs, top = self.top_module))
        verbose_append("open_mw_lib {top}_LIB".format(top = self.top_module))

        # Reading other input libraries 
        fast_lib = ' '.join(self.technology.read_libs([hammer_tech.filters.timing_db_filter], HammerTechnologyUtils.to_plain_item,
                            extra_pre_filters=[self.filter_for_lib_type("fast", "fast")]))
        slow_lib = ' '.join(self.technology.read_libs([hammer_tech.filters.timing_db_filter], HammerTechnologyUtils.to_plain_item,
                            extra_pre_filters=[self.filter_for_lib_type("slow", "slow")]))

        verbose_append("set_app_var link_library {slow_lib}".format(slow_lib = slow_lib))
        verbose_append("set_app_var target_library {slow_lib}".format(slow_lib = slow_lib))
        verbose_append("set_min_library {slow_lib} -min_version {fast_lib}".format(slow_lib = slow_lib, fast_lib = fast_lib))
        verbose_append('set symbol_library {}'.format(""))

        verbose_append("set_tlu_plus_files -max_tluplus {tlu_max} -min_tluplus {tlu_min} -tech2itf_map {tlu_map}".format(
                        tlu_max = tlu_max_caps, tlu_min = tlu_min_caps, tlu_map = tlu_map))
        
        # Reading the Gate Level Netlist and SDC files
        verbose_append("read_verilog -top {top} ../syn-rundir/results/{top}.mapped.v".format(top=self.top_module))                 
        verbose_append("read_sdc ../syn-rundir/results/{top}.mapped.sdc".format(top=self.top_module))

        # Future TODO: Implement MCMM flow
        # Read the required inputs for writing mcmm.scenarios file
        # Create a function to write mcmm.scenarios file
        # Write commands to carry out mcmm analysis on the scenarios in mcmm.scenarios file
        
        grid_unit = self.technology.get_grid_unit()
        verbose_append("set_user_grid -user_grid {o}{o}0 0{c} {o}{grid_unit} {grid_unit}{c}{c}".format(
                       o='{', c='}', grid_unit = grid_unit))

        return True

    def floorplan_design(self) -> bool:
        """ Generates the floorplan based on the floorplan mode specified in the Hammer IR """
        floorplan_tcl = os.path.join(self.run_dir, "floorplan.tcl")
        floorplan_mode = self.get_setting("par.icc.floorplan_mode")
        
        if floorplan_mode == "annotations":
            raise NotImplementedError("Not implemented")
        elif floorplan_mode == "manual":
            floorplan_script = self.get_setting("par.icc.floorplan_script")
            if floorplan_script == "null" or floorplan_script == "":
                self.logger.error("floorplan_mode is manual but no floorplan_script specified")
                return False
            copyfile(floorplan_script, floorplan_tcl)
        elif floorplan_mode == "generate":
            with open(floorplan_tcl, "w") as f:
                f.write("\n".join(self.generate_floorplan_tcl()))    
        elif floorplan_mode == "default":
            with open(floorplan_tcl, "w") as f:
                f.write("create_floorplan -control_type aspect_ratio -core_aspect_ratio 1 -core_utilization 0.7 -left_io2core 0 "
                        "-bottom_io2core 0 -right_io2core 0 -top_io2core 0 -start_first_row")                        
        else:
            self.logger.error("Invalid floorplan_mode %s" % (floorplan_mode))
            return False

        with open(floorplan_tcl, "a") as f:
            f.write("\nconnect_net {o}{power_net}{c} {o}{power_port}{c}".format(
                     o='{', power_net=self.get_setting("par.icc.MW_POWER_NET"),
                     power_port=self.get_setting("par.icc.MW_POWER_PORT"), c='}'))
            f.write("\nconnect_net {o}{ground_net}{c} {o}{ground_port}{c}".format(
                     o='{', ground_net=self.get_setting("par.icc.MW_GROUND_NET"),
                     ground_port=self.get_setting("par.icc.MW_GROUND_PORT"), c='}'))

            
        self.verbose_append("source -echo -verbose {}".format(floorplan_tcl))
        return True

    def generate_floorplan_tcl(self) -> List[str]:
        """
        Generates the Floorplan Tcl sript based on the placement constraints specified by user in Hammer IR
        """
        output = [] #type: List[str]

        global_top_layer = self.get_setting("par.blockage_spacing_top_layer")  # type: Optional[str]	
        floorplan_constraints = self.get_placement_constraints()
        	
        for constraint in floorplan_constraints:
        		
            new_path = "/".join(constraint.path.split("/")[1:])
        		
            if new_path == "":
                assert constraint.type == PlacementConstraintType.TopLevel, "Top must be a top-level/chip size constraint"
                margins = constraint.margins
                assert margins is not None
                # TCL Command for top-level chip dimensions
                output.append("create_floorplan -control_type width_and_height -core_aspect_ratio 1.0 -core_utilization 0.7 " 
                              "-core_width {width} -core_height {height} -left_io2core {left} -bottom_io2core {bottom} "  
                              "-right_io2core {right} -top_io2core {top} -start_first_row".format(width=constraint.width, 
                               height=constraint.height, left=margins.left, bottom=margins.bottom, right=margins.right, top=margins.top))
                       
            else:
                orientation = constraint.orientation if constraint.orientation is not None else "r0" # need to find if orientation is relevant 
                                                                                                     # in Synopsys IC Compiler
               
                if constraint.create_physical:
                    # TODO: need to find TCL command equivalent to "create_inst" for Synopsys IC Compiler
                    pass

                if constraint.type == PlacementConstraintType.Dummy: #does nothing with the constraint
                    pass
                
                elif constraint.type == PlacementConstraintType.Placement:
                    # TO VERIFY/TEST
                    output.append("create_bounds -name {bound_name} -coordinate {{{llx1} {lly1} {urx1} {ury1}}} -type hard "
                                  "{inst_name}".format(bound_name=new_path, llx1=constraint.x, lly1=constraint.y,
                                   urx1=constraint.x + constraint.width, ury1=constraint.y + constraint.height, inst_name=new_path))
                
                elif constraint.type in [PlacementConstraintType.HardMacro, PlacementConstraintType.Hierarchical]:
                    output.append("create_bounds -name {bound_name} -coordinate {{{llx1} {lly1} {urx1} {ury1}}} -type hard "
                                  "{inst_name}".format(bound_name=new_path, llx1=constraint.x, lly1=constraint.y,
                                   urx1=constraint.x + constraint.width, ury1=constraint.y + constraint.height, inst_name=new_path))
                    spacing = self.get_setting("par.blockage_spacing")
                    if constraint.top_layer is not None:
                        current_top_layer = constraint.top_layer
                    elif global_top_layer is not None:
                        current_top_layer = global_top_layer
                    else:
                        current_top_layer = None
                    if current_top_layer is not None:
                        bot_layer = self.get_stackup().get_metal_by_index(1).name
                        output.append("set_keepout_margin -outer {{{s} {s} {s} {s}}} {inst_name}".format(s=spacing, inst_name=new_path))
                        # ICC doesn't have an option of specifying a route halo where the cost of routing is very high
                        # in the corridor around the macros/blocks
               
                elif constraint.type == PlacementConstraintType.Obstruction:
                    # TO VERIFY/TEST
                    obs_types = get_or_else(constraint.obs_types, [])  # type: List[ObstructionType]
                    if ObstructionType.Place in obs_types:
                        output.append("create_placement_blockage -bbox {{{llx1} {lly1} {urx1} {ury1}}} -type hard".format(llx1=constraint.x,
                                       lly1=constraint.y, urx1=constraint.x + constraint.width, ury1=constraint.y + constraint.height))

                    if ObstructionType.Route in obs_types:
                        output.append("create_route_guide -coordinate {{{llx1} {lly1} {urx1} {ury1}}} -no_signal_layers {{{layers}}} " 
                                      "-zero_min_spacing".format(llx1=constraint.x, lly1=constraint.y, urx1=constraint.x + 
                                                                 constraint.width, ury1=constraint.y + constraint.height,
                                        layers="all" if constraint.layers is None else " ".join(get_or_else(constraint.layers, [])) ))

                    if ObstructionType.Power in obs_types:
                        output.append("create_route_guide -coordinate {{{llx1} {lly1} {urx1} {ury1}}} -no_preroute_layers {{{layers}}} " 
                                      "-zero_min_spacing".format(llx1=constraint.x, lly1=constraint.y, urx1=constraint.x +
                                                                 constraint.width, ury1=constraint.y + constraint.height,
                                        layers="all" if constraint.layers is None else " ".join(get_or_else(constraint.layers, [])) ))

                else:
                    assert False, "Should not reach here"
                		 					 
        return output    
   
    def place_bumps(self) -> bool:
        """ Produces Bump Placement Tcl script based on the Bumps Mode specified in the Hammer IR """
        place_bumps_tcl = os.path.join(self.run_dir, "place_bumps.tcl")

        bumps_mode = self.get_setting("vlsi.inputs.bumps_mode")
        if bumps_mode == "manual":
            bumps = self.get_bumps()
            if bumps is not None:  
                # Call the Bump-Placement API method
                with open(place_bumps_tcl, "w") as f:
                        f.write("\n".join(self.generate_bumps_tcl()))
                self.verbose_append("source -echo -verbose {}".format(place_bumps_tcl))
            else:
                self.logger.error("Bumps Placement Mode is manual but no bumps are specified!!")
                return False
        elif bumps_mode == "empty":
            self.logger.warning("Bumps Placement Mode is empty. No bumps are added.")
        else:
            self.logger.error("Invalid Bumps Placement Mode %s" % (bumps_mode))
            return False

        return True   

    def generate_bumps_tcl(self) -> List[str]:
        """
        Generates the Bump Placement Tcl sript based on the bumps assignment constraints specified by the user in Hammer IR
        """
        output = [] #type: List[str]
       
        bumps = self.get_bumps()
        assert isinstance(bumps, BumpsDefinition) # bumps not be None
        bump_array_width = Decimal(str((bumps.x - 1) * bumps.pitch))
        bump_array_height = Decimal(str((bumps.y - 1) * bumps.pitch))
        fp_consts = self.get_placement_constraints()
        fp_width = Decimal(0)
        fp_height = Decimal(0)
           
        for const in fp_consts:
            if const.type == PlacementConstraintType.TopLevel:
                fp_width = const.width
                fp_height = const.height
        if fp_width == 0 or fp_height == 0:
            raise ValueError("Floorplan does not specify a TopLevel constraint or it has no dimensions")
           
        # Center bump array in the middle of floorplan
        bump_offset_x = (Decimal(str(fp_width)) - bump_array_width) / 2
        bump_offset_y = (Decimal(str(fp_height)) - bump_array_height) / 2
        power_ground_nets = list(map(lambda x: x.name, self.get_independent_power_nets() + self.get_independent_ground_nets()))        
        block_layer = self.get_setting("vlsi.technology.bump_block_cut_layer")  # type: str
            
        for bump in bumps.assignments:
            output.append("place_flip_chip_array -physical_lib_cell {bump_cell} -prefix \"bump_{c}.{r}\" -start_point {{{x} {y}}} "
                          "-number 1 -delta {{{pitch} {pitch}}} -repeat {{1 1}} -cell_origin center".format(
                          bump_cell=bump.custom_cell if bump.custom_cell is not None else bumps.cell, 
                          c=bump.x, r=bump.y, x=bump_offset_x + Decimal(str(bump.x - 1)) * Decimal(str(bumps.pitch)),
                          y=bump_offset_y + Decimal(str(bump.y - 1)) * Decimal(str(bumps.pitch)), pitch=bumps.pitch))
                
            if not bump.no_connect:
                if bump.name in power_ground_nets:
                    output.append("change_selection [get_cells -all bump_{c}.{r}*]".format(c=bump.x, r=bump.y))
                    output.append('set_flip_chip_type -personality "PG" [get_selection]')
                    output.append('change_selection ""')
                else:
                    output.append("change_selection [get_cells -all bump_{c}.{r}*]".format(c=bump.x, r=bump.y))
                    output.append('set_flip_chip_type -personality "Signal" [get_selection]')
                    output.append('change_selection ""')
                     
            output.append("create_route_guide -coordinate {{{llx1} {lly1} {urx1} {ury1}}} -no_signal_layers {{{layer}}} "
                          "-zero_min_spacing".format(llx1="get_attribute [get_cells -all bump_{c}.{r}*] bbox_llx".format(c=bump.x, r=bump.y),                                           lly1="get_attribute [get_cells -all bump_{c}.{r}*] bbox_lly".format(c=bump.x, r=bump.y), 
                                                     urx1="get_attribute [get_cells -all bump_{c}.{r}*] bbox_urx".format(c=bump.x, r=bump.y), 
                                                     ury1="get_attribute [get_cells -all bump_{c}.{r}*] bbox_ury".format(c=bump.x, r=bump.y), 
                                                     layer=block_layer))
            
        # Assigning the created bumps to nets
        ## TODO: To find a way to incorporate variable no. of PG and Signal drivers
        ## TODO: To remove the dependency on the below ref supply dependent keys
        vdd_ref = self.get_setting("par.icc.VDD_ref") 
        vddio_ref = self.get_setting("par.icc.VDDIO_ref")
        vss_ref = self.get_setting("par.icc.VSS_ref")
        vssio_ref = self.get_setting("par.icc.VSSIO_ref")
        signal_ref = self.get_setting("par.icc.Signal_ref")
    
        # Select the P/G drivers based on the ref cell name
        output.append('change_selection [get_cells -all -hierarchical -filter {{ref_name="{VDD}" || ref_name="{VDDIO}" || ref_name="{VSS}" || ref_name="{VSSIO}"}}]'.format(
            VDD=vdd_ref, VDDIO=vddio_ref, VSS=vss_ref, VSSIO=vssio_ref))
        output.append('set_flip_chip_type -personality "PG" [get_selection]') 
            
        # Select the Signal drivers based on the ref cell name
        output.append('change_selection [get_cells -all -hierarchical -filter {{ref_name="{Signal}"}}]'.format(Signal=signal_ref))
        output.append('set_flip_chip_type -personality "Signal" [get_selection]')
            
        output.append("assign_flip_chip_nets")
            
        return output
    
    def place_tap_cells(self) -> bool:
        """
        Generates the ICC Tcl sript for placement of tapcells based on the tapcell constraints specified by the user in Hammer IR
        """
        tap_cells = self.technology.get_special_cell_by_type(CellType.TapCell)

        if len(tap_cells) == 0:
            self.logger.warning("Tap cells are improperly defined in the tech plugin and will not be added."
                                "This step should be overridden with a user hook.")
            return True

        tap_cell = tap_cells[0].name[0]

        try:
            interval = self.get_setting("vlsi.technology.tap_cell_interval")
            offset = self.get_setting("vlsi.technology.tap_cell_offset")
            self.append("add_tap_cell_array -master_cell_name {tap_cell} -distance {cell_interval} -offset {row_offset} "
                        "-pattern normal".format(tap_cell=tap_cell, cell_interval=interval, row_offset=offset))
        except KeyError:
            pass
        finally:
            self.logger.warning("You have not overridden place_tap_cells. By default this step adds a simple set of tapcells or "
                                "does nothing; you will have trouble with power strap creation later.")
        return True
     
    def place_pins(self) -> bool:
        """ Generates Pin Placement Tcl script based on the Pin Mode specified in Hammer IR """
        place_pins_tcl = os.path.join(self.run_dir, "place_pins.tcl")
        
        pin_mode = self.get_setting("vlsi.inputs.pin_mode")
        if pin_mode == "none":
            # This mode lets the CAD tool do whatever it wants (providing no constraints).
            # Typically this is sane but unpredictable.
            self.logger.warning("Pin Placement Mode is none!! CAD Tool does whatever it wants to.")
        elif pin_mode == "generated":             
            pin_generation_mode = self.get_setting("vlsi.inputs.pin.generate_mode")
            if (pin_generation_mode == "full_auto" or pin_generation_mode == "semi_auto"):    
                pins = self.get_pin_assignments()
                if pins is not None:
                    # Call the Pin-Placement API method
                    with open(place_pins_tcl, "w") as f:
                            f.write("\n".join(self.generate_pins_tcl()))
                    self.verbose_append("source -echo -verbose {}".format(place_pins_tcl))
                else:
                    self.logger.error("Pins are set to be generated but no pins are specified!!!")
                    return False
            else:
                self.logger.error("Invalid Pin Generation Mode %s" % (pin_generation_mode))
                return False
        elif pin_mode == "auto":
            # (not implemented yet) looks at the hierarchy and guess where pins should go
            raise NotImplementedError("Automatic Pin Placement Mode not implemented yet")
        else:
            self.logger.error("Invalid Pin Placement Mode %s" % (pin_mode))
            return False

        return True    

    def generate_pins_tcl(self) -> List[str]:
        """Generates the Pin Placement Tcl sript based on the pin assignment constraints specified by the user in Hammer IR """
        output = [] #type: List[str]

        output.append("set_fp_pin_constraints -block_level -hard_constraints {layer location} -use_physical_constraints on")
        output.append("create_port {o}{power_port} {ground_port}{c} -direction inout".format(
                     o='{', power_port=self.get_setting("par.icc.MW_POWER_PORT"),
                     ground_port=self.get_setting("par.icc.MW_GROUND_PORT"), c='}'))

        fp_consts = self.get_placement_constraints()
        topconst = None  # type: Optional[PlacementConstraint]
        for const in fp_consts:
            if const.type == PlacementConstraintType.TopLevel:
                topconst = const
        assert topconst is not None, "Cannot find top-level constraints to place pins"

        const = cast(PlacementConstraint, topconst)
        assert isinstance(const.margins, Margins), "Margins must be defined for the top level"
        fp_llx = const.margins.left
        fp_lly = const.margins.bottom
        fp_urx = const.width - const.margins.right
        fp_ury = const.height - const.margins.top

        pin_assignments = self.get_pin_assignments() 

        for pin in pin_assignments:
                # Editing the features of the pin : width, depth, location, side, layer
                synopsys_side_num = None
                side_arg = ""
                if pin.side is not None:
                    if pin.side == "internal":
                        side_arg = "-off_edge location"
                    else:
                        side_map = {"left": 1, "top": 2, "right": 3, "bottom": 4}
                        synopsys_side_num = side_map[str(pin.side)]
                        side_arg = "-side {side_num}".format(side_num=synopsys_side_num)                        
 
                location_arg = ""
                if pin.location is None:
                    pass # start and end points of pin placement along edges is automatically handled by Synopsys ICC
                else:
                   location_arg = "-location {{ {x} {y} }}".format(x=pin.location[0], y=pin.location[1])

                layers_arg = ""
                if pin.layers is not None and len(pin.layers) > 0: # pin.layers must consist of consecutive metal layers only 
                    layers_arg = "-layers {{ {} }}".format(" ".join(pin.layers))

                width_arg = get_or_else(optional_map(pin.width, lambda f: "-width {f}".format(f=f)), "")
                depth_arg = get_or_else(optional_map(pin.depth, lambda f: "-depth {f}".format(f=f)), "")

                cmd = [
                       "set_pin_physical_constraints",
                       "-pin_name", pin.pins,
                       side_arg,
                       location_arg,
                       layers_arg,
                       width_arg,
                       depth_arg
                      ]

                output.append(" ".join(cmd))

        return output

    def power_straps(self) -> bool:
        """Place the power straps for the design."""
        power_straps_tcl = os.path.join(self.run_dir, "power_straps.tcl")
        power_straps_mode = self.get_setting("par.power_straps_mode")
        
        if power_straps_mode == "empty":
            self.logger.warning("Power Straps Mode is empty!! No power straps are generated")
        elif power_straps_mode == "manual":
            # Reads the manual power straps script contents and copies it into power_straps_tcl
            power_straps_script = self.get_setting("par.power_straps_script_contents")
            copyfile(power_straps_script, power_straps_tcl)     
        elif power_straps_mode == "generate":
            # Generate the Power Straps by passing the input from Hammer IR to the Power Straps Generation API
            # Call the Power Straps API methods
            with open(power_straps_tcl, "w") as f:
                    f.write("\n".join(self.create_power_straps_tcl()))
        else:
            self.logger.error("Invalid Power Straps Mode %s" % (power_straps_mode))
            return False

        self.verbose_append("source -echo -verbose {}".format(power_straps_tcl))
        return True

    def specify_std_cell_power_straps(self, blockage_spacing: Decimal, bbox: Optional[List[Decimal]], nets: List[str]) -> List[str]:
        # TODO: resolve doubts regarding via
        """
        Generate a list of TCL commands that build the low-level standard cell power strap rails.
        The layer is set by technology.core.std_cell_rail_layer, which should be the highest metal layer in the std cell rails.
       
        :param bbox: The optional (2N)-point bounding box of the area to generate straps. By default the entire core area is used.
        :param nets: A list of power net names (e.g. ["VDD", "VSS"]). Currently only two are supported.
        :return: A list of TCL commands that will generate power straps on rails.        

        """
        layer_name = self.get_setting("technology.core.std_cell_rail_layer")
        layer = self.get_stackup().get_metal(layer_name)
        master = self.get_setting("technology.core.tap_cell_rail_reference") 

        results = [
            "# Power strap definition for layer {} (rails):\n".format(layer_name),
            "set_preroute_special_rules -name splrule -pin_layer {layer}"
            " -macro_cells {{{cells}}}".format(layer=layer_name, cells=master)
            # need to find how to set bottom and top via stack layers in Synopsys ICC    
                  ]
         
        options = [
            "-layer", layer_name,
            "-configure", "macros",
            "-special_rules", "splrule",
            "-direction", str(layer.direction),
            "-nets", "{ %s }" % " ".join(nets)
        ]
       
        ## Note: The area inside which power straps will be generated can't be controlled in Synopsys ICC. The workaround is to use
        ##       parameters like -start and -stop of the create_power_straps command. By default, the entire core area will be used. 
       
        results.append("create_power_straps " + " ".join(options) + "\n")
        return results

    def specify_power_straps(self, layer_name: str, bottom_via_layer_name: str, blockage_spacing: Decimal, pitch: Decimal, width: Decimal, spacing: Decimal, offset: Decimal, bbox: Optional[List[Decimal]], nets: List[str], add_pins: bool) -> List[str]:
        # TODO: resolve doubts regarding via
        """
        Generate a list of TCL commands that will create power straps on a given layer.
        This is a low-level, cad-tool-specific API. It is designed to be called by higher-level methods, so calling this directly is not
        recommended.
        This method assumes that power straps are built bottom-up, starting with standard cell rails.
       
        :param layer_name: The layer name of the metal on which to create straps.
        :param bottom_via_layer_name: The layer name of the lowest metal layer down to which to drop vias.
        :param blockage_spacing: The minimum spacing between the end of a strap and the beginning of a macro or blockage.
        :param pitch: The pitch between groups of power straps (i.e. from left edge of strap A to the next left edge of strap A).
        :param width: The width of each strap in a group.
        :param spacing: The spacing between straps in a group.
        :param offset: The offset to start the first group.
        :param bbox: The optional (2N)-point bounding box of the area to generate straps. By default the entire core area is used.
        :param nets: A list of power nets to create (e.g. ["VDD", "VSS"], ["VDDA", "VSS", "VDDB"],  ... etc.).
        :param add_pins: True if pins are desired on this layer; False otherwise.
        :return: A list of TCL commands that will generate power straps.        
  
        """
        results = ["# Power strap definition for layer %s:\n" % layer_name]
        # results.extend([ need to find how to set bottom and top via stack layers, antenna rule and blockage spacing ])
        
        layer = self.get_stackup().get_metal(layer_name)

        fp_consts = self.get_placement_constraints()
        fp_width = Decimal(0)
        fp_height = Decimal(0)
           
        for const in fp_consts:
            if const.type == PlacementConstraintType.TopLevel:
                fp_width = const.width
                fp_height = const.height
        if fp_width == 0 or fp_height == 0:
            raise ValueError("Floorplan does not specify a TopLevel constraint or it has no dimensions")
       
        options = [
            "-direction", str(layer.direction),
            "-layer", layer_name,
            "-nets", "{%s}" % " ".join(nets),
            "-configure", "step_and_stop",
            "-step", str(pitch),
            "-pitch_within_group", str(spacing),
            "-stop", str(fp_width) if layer.direction == RoutingDirection.Vertical else str(fp_height),
            "-width", str(width),
            "-start_at", str(offset)
        ]

        # This is temporary implementation - TODO: need to find how to make a stripe into a LEF pin shape
        if (add_pins):
            options.extend([
                 "-extend_low_ends", "to_boundary_and_generate_pins",
                 "-extend_high_ends", "to_boundary_and_generate_pins"
                          ])
        else:
            pass

        index = 0
        if layer.direction == RoutingDirection.Horizontal:
            index = 1
        elif layer.direction != RoutingDirection.Vertical:
            raise ValueError("Cannot handle routing direction {d} for layer {l} when creating power straps".format(d=str(layer.direction),
                                                                                                                   l=layer_name))

        results.append("create_power_straps " + " ".join(options) + "\n")
        return results            
   
    def place_opt_design(self) -> bool:
        """Place the design and do pre-routing optimization."""
        verbose_append = self.verbose_append
        verbose_append("set_auto_disable_drc_nets -constant false")
        verbose_append("set_fix_hold [get_clocks *]")
        verbose_append("set_app_var enable_recovery_removal_arcs true")
        verbose_append("set_place_opt_strategy -layer_optimization_effort high")
        verbose_append("set_place_opt_strategy -consider_routing true")
        verbose_append("create_fp_placement -effort high")
        verbose_append("\n".join(self.pg_connection()))
        verbose_append("place_opt -effort {effort} -congestion".format(effort = self.get_setting("par.icc.place_opt_effort")))
        return True

    def clock_tree(self) -> bool:
        """Setup and route a clock tree for clock nets."""
        verbose_append = self.verbose_append
        verbose_append("check_legality -verbose")
        verbose_append("refine_placement -num_cpus 0")
        verbose_append("legalize_placement -effort high")
        verbose_append("report_placement_utilization")
        verbose_append("set_auto_disable_drc_nets -constant false")
        verbose_append("remove_ideal_network -all")
        verbose_append("set_propagated_clock [all_clocks]")
        verbose_append("set_fix_hold [get_clocks *]")
        verbose_append("set_route_zrt_common_options -read_user_metal_blockage_layer true")
        verbose_append("set_app_var timing_enable_multiple_clocks_per_reg true")
        verbose_append("set_app_var cto_enable_drc_fixing true")
        verbose_append("set cts_instance_name_prefix cts")
        verbose_append("set_si_options -delta_delay true -route_xtalk_prevention true -static_noise true")
        verbose_append("\n".join(self.pg_connection()))
        verbose_append("clock_opt -fix_hold_all_clocks -operating_condition min_max -congestion")
        return True

    def cts_opt(self) -> bool:
        """ Post CTS Optimization Script"""
        verbose_append = self.verbose_append
        verbose_append("remove_ideal_network [get_clocks *]")
        verbose_append("set_propagated_clock [all_clocks]")
        verbose_append("set_fix_hold [get_clocks *]")
        verbose_append("psynopt -only_hold_time")
        return True

    def add_fillers(self) -> bool:
        """Add filler cells"""
        stdfillers = self.technology.get_special_cell_by_type(CellType.StdFiller)

        if len(stdfillers) == 0:
            self.logger.warning(
                "The technology plugin 'special cells: stdfiller' field does not exist. It should specify a list of (non IO) filler cells. " 
                "No filler will be added. You can override this with an add_fillers hook if you do not specify filler cells in the "
                "technology plugin.")
        else:
            stdfiller = stdfillers[0].name
            filler_str = ""
            for cell in stdfiller:
                filler_str += str(cell) + ' '
            self.append("insert_stdcell_filler -cell_without_metal \"{FILLER}\" -cell_with_metal \"{FILLER}\" "
                        "-connect_to_power {{ {VDD} }} -connect_to_ground {{ {VSS} }}".format(FILLER = filler_str,
                        VDD = self.get_setting("par.icc.MW_POWER_PORT"), VSS = self.get_setting("par.icc.MW_GROUND_PORT")))    
        return True   

    def route_design(self) -> bool:
        """Routing and Post-route Optimization of the design."""
        verbose_append = self.verbose_append
        verbose_append("set_app_var enable_recovery_removal_arcs true")
        verbose_append("set_fix_hold [get_clocks *]")
        verbose_append("set_si_options -delta_delay true -static_noise true")
        verbose_append("set_route_mode_options -zroute true")
        verbose_append("set_route_zrt_global_options -default true")
        verbose_append("set_route_zrt_global_options -effort high")
        verbose_append("set_route_zrt_track_options -default true")
        verbose_append("set_app_var routeopt_drc_over_timing true")
        verbose_append("set_route_zrt_detail_options -timing_driven false -drc_convergence_effort_level high")
        verbose_append("set_route_opt_strategy -fix_hold_mode all")
        verbose_append("route_opt -effort high")
        verbose_append("route_opt -only_hold_time -xtalk_reduction -effort high -incremental")
        return True

    def write_netlist(self) -> bool:
        """
        Writing Verilog outputs after PAR
        """
        self.verbose_append("change_names -rules verilog -hierarchy")

        self.verbose_append("write_verilog -pg -no_physical_only_cells -no_corner_pad_cells -no_pad_filler_cells -no_core_filler_cells "                                  "-supply_statement all {netlist}".format(netlist = self.output_netlist_filename))

        self.verbose_append("write_verilog -no_physical_only_cells -no_pg_pin_only_cells -no_corner_pad_cells {netlist}".format(
                            netlist = self.output_sim_netlist_filename))                             

        return True

    def write_gds(self) -> bool:
        """ Write GDS after PAR """

        map_file = get_or_else(
            optional_map(self.get_gds_map_file(), lambda f: "-map_layer {}".format(f)),
            ""
        )
        
        if (self.get_setting("par.inputs.gds_merge") == "true") or (self.get_setting("par.inputs.gds_precision_mode") == "manual"):
            self.logger.warning("Synopsys ICC doesn't support gds_merge and gds_precision features."
                                "They have to be done in the later stages during Physical Verification.")

        self.verbose_append("set_write_stream_options -skip_ref_lib_cells")
        self.verbose_append("set_write_stream_options -output_pin text -keep_data_type")
        self.verbose_append("set_write_stream_options -output_pin geometry")
        self.verbose_append("set_write_stream_options {gds_map}".format(gds_map = map_file))
        self.verbose_append("write_stream -format gds -cells {top} {gds_file}".format(
                             top=self.top_module, gds_file=self.output_gds_filename))

        return True

    def write_sdf(self) -> bool:
        """
         Output the Standard Delay Format File for use in timing annotated gate level simulations
        """
        self.verbose_append("write_sdf {sdf_file}".format(sdf_file = self.output_sdf_path))
        return True

    def write_def(self) -> bool:
        """
         Output the Design Exchange Format file for use across various tools
        """ 
        self.verbose_append("write_def -output {def_file}".format(def_file = self.output_def_path))
        return True   

    def write_spefs(self) -> bool:
        """
        Output a SPEF file that contains the parasitic extraction results
        """
        self.verbose_append("extract_rc -coupling_cap")
        # For MMMC designs, ICC automatically uses the name of the TLUPlus file and the temperature associated with the corner,
        # along with the specified file name, to derive the file name of the parasitic files. 
        # The naming convention used by the tool is <tluplus_file_name>_<temperature>.<output_file_name>
        self.verbose_append("write_parasitics -format SPEF -output {run_dir}/{top}.par.spef".format(
                            run_dir=self.run_dir, top=self.top_module))
        return True

    def pg_connection(self) -> List[str]:
        """ Connects the P/G pins to the P/G nets. Synopsys ICC requires this to be executed everytime a new element is added to the layout, 
            preferably after every step of the PAR flow. """
        output = [] # List[str]
        output.append("derive_pg_connection -power_net {o}{p_net}{c} -ground_net {o}{g_net}{c}" 
                      " -power_pin {o}{p_pin}{c} -ground_pin {o}{g_pin}{c}".format(o='{',
                       p_net = self.get_setting("par.icc.MW_POWER_NET"), c='}', g_net = self.get_setting("par.icc.MW_GROUND_NET"),
                       p_pin = self.get_setting("par.icc.MW_POWER_PORT"), g_pin = self.get_setting("par.icc.MW_GROUND_PORT") ))
        output.append("derive_pg_connection -power_net {o}{p_net}{c} -ground_net {o}{g_net}{c} -tie".format(o='{',
                       p_net = self.get_setting("par.icc.MW_POWER_NET"), c='}', g_net = self.get_setting("par.icc.MW_GROUND_NET") ))
        return output                              

    @property
    def output_innovus_lib_name(self) -> str:
        return "{top}_FINAL".format(top=self.top_module)

    @property
    def generated_scripts_dir(self) -> str:
        return os.path.join(self.run_dir, "generated-scripts")

    @property
    def open_chip_script(self) -> str:
        return os.path.join(self.generated_scripts_dir, "open_chip")

    @property
    def open_chip_tcl(self) -> str:
        return self.open_chip_script + ".tcl"

    # Future TODO: To include write_ilm once hierarchical flow is implemented

    def write_design(self) -> bool:
        # Save the ICC design.
        self.verbose_append("save_mw_cel -as {lib_name}".format(lib_name=self.output_innovus_lib_name))

        # Write netlist
        self.write_netlist()

        # GDS streamout.
        self.write_gds()

        # Write DEF
        self.write_def()
 
        # Write SDF
        self.write_sdf()
        

        # Write SPEF
        self.write_spefs()

        # Make sure that generated-scripts exists.
        os.makedirs(self.generated_scripts_dir, exist_ok=True)

        self.ran_write_design=True

        return True

    def run_icc(self) -> bool:
        # Quit ICC
        self.verbose_append("exit")

        # Create par script.
        par_tcl_filename = os.path.join(self.run_dir, "par.tcl")
        icc_log_file = os.path.join(self.run_dir, "icc.log")
        with open(par_tcl_filename, "w") as f:
            f.write("\n".join(self.output))

        # Make sure that generated-scripts exists.
        os.makedirs(self.generated_scripts_dir, exist_ok=True)

        # Create open_chip script pointing to latest (symlinked to post_<last ran step>).
        with open(self.open_chip_tcl, "w") as f:
            f.write("""
        cd {run_dir}
        open_mw_lib -readonly {top}_LIB
        open_mw_cel -readonly latest
                """.format(run_dir=self.run_dir, top=self.top_module))
        
        with open(self.open_chip_script, "w") as f:
            f.write("""#!/bin/bash
        cd {run_dir}
        source enter
        $ICC_BIN -gui -f {open_chip_tcl}
                """.format(run_dir=self.run_dir, open_chip_tcl=self.open_chip_tcl))
        os.chmod(self.open_chip_script, 0o755)
         
        # Build args.
        args = [
            self.get_setting("par.icc.icc_bin"),
            "-no_gui",  # Prevent the GUI popping up.
            "-f", par_tcl_filename,
            "-output_log_file", icc_log_file
        ]

        # Temporarily disable colours/tag to make run output more readable.
        # TODO: think of a more elegant way to do this?
        HammerVLSILogging.enable_colour = False
        HammerVLSILogging.enable_tag = False
        self.run_executable(args, cwd=self.run_dir)  # TODO: check for errors and deal with them
        HammerVLSILogging.enable_colour = True
        HammerVLSILogging.enable_tag = True

        # TODO: check that par run was successful

        return True
    
    def filter_for_lib_type(self, nmos, pmos) -> Callable[[hammer_tech.Library],bool]:
        """
        Selecting libraries that match given nmos and pmos type
        """
        def extraction_func(lib: hammer_tech.Library) -> bool:
            if lib.corner is None or lib.corner.nmos is None or lib.corner.pmos is None:
                return False
            nmos_type = str(lib.corner.nmos)
            pmos_type = str(lib.corner.pmos)
            if nmos_type == nmos:
                if pmos_type == pmos:
                    return True
                else:
                    return False
            else:
                return False
        return extraction_func   

def icc_global_settings(ht: HammerTool) -> bool:
    """Settings that need to be reapplied at every tool invocation"""
    assert isinstance(ht, HammerPlaceAndRouteTool)
    assert isinstance(ht, SynopsysTool)
    ht.create_enter_script()

    # Python sucks here for verbosity
    verbose_append = ht.verbose_append

    # Generic settings
    # Additional generic settings needed at every tool invocation can be added if needed
    verbose_append("set_host_options -max_cores {}".format(ht.get_setting("vlsi.core.max_threads")))

    return True

tool = ICC
