#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  hammer-vlsi plugin for Synopsys ICC.
#
#  See LICENSE for licence details.

import os
from shutil import copyfile
from os.path import dirname
from typing import List, Optional,Callable, Tuple, Set, Any, cast, Dict
from decimal import Decimal

import hammer_tech
from hammer_tech import RoutingDirection, Metal
from hammer_vlsi import HammerPlaceAndRouteTool, HammerToolStep, PlacementConstraintType, HierarchicalMode, ObstructionType, Margins, Supply, PlacementConstraint, MMMCCornerType
from hammer_vlsi import SynopsysTool
from hammer_logging import HammerVLSILogging
from hammer_tech import HammerTechnologyUtils
import specialcells
from specialcells import CellType, SpecialCell

class ICC(HammerPlaceAndRouteTool, SynopsysTool):
    def fill_outputs(self) -> bool:
        # TODO: implement
        return True
    
    def place_pins(self) -> List[str]:
        """
        Generates the Pin Placement Tcl sript based on the pin assignment constraints specified by the user in Hammer IR
        """
        output = [] #type: List[str]
       
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

        promoted_pins = []  # type: List[str]
        for pin in pin_assignments:
            if pin.preplaced:
                ## TODO First set promoted pins to handle preplaced pins feature of the API
                # find equivalent for set_promoted_macro_pin
                pass
            else:
                # Editing the features of the pin : width, depth, location, side, layer
                
                synopsys_side_num = None
                side_arg = ""
                if pin.side is not None:
                    if pin.side == "internal":
                        side_arg = "-off_edge location"
                    else:
                        side_map = {"left": 1, "top": 2, "right": 3, "bottom": 4}
                        synopsys_side_num = side_map(str(pin.side))
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
                output.append("set_fp_pin_constraints -block_level -hard_constraints {layer location} -use_physical_constraints on")
                output.append("place_fp_pins -block_level")
                
        for pin in promoted_pins:
            # find equivalent for assign_io_pins
            pass

        return output

    
    def place_bumps(self) -> List[str]:
        """
        Generates the Bump Placement Tcl sript based on the bumps assignment constraints specified by the user in Hammer IR
        """
        output = [] #type: List[str]
       
        bumps = self.get_bumps()
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
            
        pg_bumps = [] # type: List[str]
        signal_bumps = [] # type: List[str]
        for bump in bumps.assignments:
            output.append("place_flip_chip_array -physical_lib_cell {bump_cell} -prefix \"bump_{c}.{r}\" -start_point {{x} {y}} "
                          "-number 1 -delta {{pitch} {pitch}} -repeat {1 1} -cell_origin center".format(
                          bump_cell=bump.custom_cell if bump.custom_cell is not None else bumps.cell, 
                          c=bump.x, r=bump.y, x=bump_offset_x + Decimal(str(bump.x - 1)) * Decimal(str(bumps.pitch)),
                          y=bump_offset_y + Decimal(str(bump.y - 1)) * Decimal(str(bumps.pitch)), pitch=bumps.pitch))
                
            if not bump.no_connect:
                if bump.name in power_ground_nets:
                    pg_bumps.append(bump.name)
                else:
                    signal_bumps.append(bump.name)
                     
            output.append("create_route_guide -coordinate {{llx1} {lly1} {urx1} {ury1}} -no_signal_layers {{{layer}}} "
                          "-zero_min_spacing".format(llx1="get_attribute [get_cells bump_{c}.{r}] bbox_llx".format(c=bump.x,r=bump.y),
                                                     lly1="get_attribute [get_cells bump_{c}.{r}] bbox_lly".format(c=bump.x, r=bump.y), 
                                                     urx1="get_attribute [get_cells bump_{c}.{r}] bbox_urx".format(c=bump.x, r=bump.y), 
                                                     ury1="get_attribute [get_cells bump_{c}.{r}] bbox_ury".format(c=bump.x, r=bump.y), 
                                                     layer=block_layer))
            
        # Assigning the created bumps to nets
        ## TODO: To find a way to incorporate variable no. of PG and Signal drivers
        vdd_ref = self.get_setting("par.icc.VDD_ref") 
        vddio_ref = self.get_setting("par.icc.VDDIO_ref")
        vss_ref = self.get_setting("par.icc.VSS_ref")
        vssio_ref = self.get_setting("par.icc.VSSIO_ref")
        signal_ref = self.get_setting("par.icc.Signal_ref")
            
        # Select the P/G drivers based on the ref cell name
        output.append('change_selection [get_cells -all -hierarchical -filter {ref_name="{VDD}" || ref_name="{VDDIO}" || ref_name="{VSS}" || ref_name="{VSSIO}"}]'.format(VDD=vdd_ref, VDDIO=vddio_ref, VSS=vss_ref, VSSIO=vssioref))
        # Select the P/G bumps
        for bump in pg_bumps:
            output.append("change_selection -add [get_cells {bump_cell}]".format(bump_cell=bump))
        output.append('set_flip_chip_type -personality "PG" [get_selection]') 
            
        # Select the Signal drivers based on the ref cell name
        output.append('change_selection [get_cells -all -hierarchical -filter {ref_name="{Signal}"}]'.format(Signal=signal_ref))
        for bump in signal_bumps:
            output.append("change_selection -add [get_cells {bump_cell}]".format(bump_cell=bump))
        output.append('set_flip_chip_type -personality "Signal" [get_selection]')
            
        output.append("assign_flip_chip_nets")
            
        return output

    def place_tap_cells(self) -> List[str]:
        """
        Generates the ICC Tcl sript for placement of tapcells based on the tapcell constraints specified by the user in Hammer IR
        """
        output = [] #type: List[str]
        tap_cells = self.technology.get_special_cell_by_type(CellType.TapCell)

        tap_cell = tap_cells[0].name[0]

        try:
            interval = self.get_setting("vlsi.technology.tap_cell_interval")
            offset = self.get_setting("vlsi.technology.tap_cell_offset")
            output.append("add_tap_cell_array -master_cell_name {tap_cell} -distance {cell_interval} -offset {row_offset} "
                          "-pattern normal".format(tap_cell=tap_cell, cell_interval=interval, row_offset=offset))
        except KeyError:
            pass
        finally:
            self.logger.warning("You have not overridden place_tap_cells. By default this step adds a simple set of tapcells or "
                                "does nothing; you will have trouble with power strap creation later.")
        return output


    def specify_power_straps(self, layer_name: str, bottom_via_layer_name: str, blockage_spacing: Decimal, pitch: Decimal, width: Decimal, spacing: Decimal, offset: Decimal, bbox: Optional[List[Decimal]], nets: List[str], add_pins: bool) -> List[str]:
        # TODO: resolve doubts regarding via and power ring settings
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
        num_strap_sets = self.get_setting("par.icc.power_straps_num_groups")
       
        options = [
            "-direction", str(layer.direction),
            "-layer", layer_name,
            "-nets", "{%s}" % " ".join(nets),
            "-configure", "groups_and_step",
            "-step", str(pitch),
            "-pitch_within_group", str(spacing),
            "-num_groups", num_strap_sets,
            "-width", str(width),
            "-start_at", str(offset)
        ]

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

    def specify_std_cell_power_straps(self, blockage_spacing: Decimal, bbox: Optional[List[Decimal]], nets: List[str]) -> List[str]:
        # TODO: resolve doubts regarding via and power ring settings
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
            " -macro_cells {{cells}}".format(layer=layer_name, cells=master)
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
       
        # if bbox is not None:
        #    options.extend([
        #        "-area", "{ %s }" % " ".join(map(str, bbox))
        #    ])
       
        results.append("create_power_straps " + " ".join(options) + "\n")
        return results

    def generate_floorplan_tcl(self) -> List[str]:
        """
        Generates the Floorplan Tcl sript based on the placement constraints specified by user in Hammer IR
        """
        output = [] #type: List[str]
        	
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
                    # TODO
                    # need to find TCL command equivalent to "create_inst" for Synopsys IC Compiler
                    pass

                if constraint.type == PlacementConstraintType.Dummy: #does nothing with the constraint
                    pass
                
                elif constraint.type == PlacementConstraintType.Placement:
                    # TO VERIFY/TEST
                    output.append("define_user_attribute -type boolean -class cell attribute_{name}".format(name=new_path))
                    output.append("set_attribute {inst_name} attribute_{name} true".format(inst_name=new_path, name=new_path))
                    output.append("create_placement_blockage -bbox {{llx1} {lly1} {urx1} {ury1}} -type partial "
                                  "-blocked_percentage 50 -category attribute_{name}".format(llx1=constraint.x, lly1=constraint.y,
                                   urx1=constraint.x + constraint.width, urx2=constraint.y + constraint.height, name=new_path))
                
                elif constraint.type in [PlacementConstraintType.HardMacro, PlacementConstraintType.Hierarchical]:
                    # TODO
                    # need to find TCL command equivalent to "place_inst" for Synopsys IC Compiler
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
                        # need to find TCL command equivalent to "create_route_halo" for Synopsys IC Compiler
               
                elif constraint.type == PlacementConstraintType.Obstruction:
                    # TO VERIFY/TEST
                    obs_types = get_or_else(constraint.obs_types, [])  # type: List[ObstructionType]
                    if ObstructionType.Place in obs_types:
                        output.append("create_placement_blockage -bbox {{llx1} {lly1} {urx1} {ury1}} -type hard".format(llx1=constraint.x,
                                       lly1=constraint.y, urx1=constraint.x + constraint.width, urx2=constraint.y + constraint.height))

                    if ObstructionType.Route in obs_types:
                        output.append("create_route_guide -coordinate {{llx1} {lly1} {urx1} {ury1}} -no_signal_layers {{{layers}}} " 
                                      "-zero_min_spacing".format(llx1=constraint.x, lly1=constraint.y, urx1=constraint.x + 
                                                                 constraint.width, urx2=constraint.y + constraint.height,
                                        layers="all" if constraint.layers is None else " ".join(get_or_else(constraint.layers, [])) ))

                    if ObstructionType.Power in obs_types:
                        output.append("create_route_guide -coordinate {{llx1} {lly1} {urx1} {ury1}} -no_preroute_layers {{{layers}}} " 
                                      "-zero_min_spacing".format(llx1=constraint.x, lly1=constraint.y, urx1=constraint.x +
                                                                 constraint.width, urx2=constraint.y + constraint.height,
                                        layers="all" if constraint.layers is None else " ".join(get_or_else(constraint.layers, [])) ))

                else:
                    assert False, "Should not reach here"
                		 					 
        return output 
  
    @property
    def steps(self) -> List[HammerToolStep]:
        return self.make_steps_from_methods([
            self.main_step
        ])

    def tool_config_prefix(self) -> str:
        return "par.icc"

    @property
    def env_vars(self) -> Dict[str, str]:
        v = dict(super().env_vars)
        v["ICC_BIN"] = self.get_setting("par.icc.icc_bin")
        return v           
 
    def main_step(self) -> bool:
        # Locate reference methodology tarball.
        synopsys_rm_tarball = self.get_synopsys_rm_tarball("ICC")
        
        icc_bin = self.env_vars["ICC_BIN"]       
 
        # Generate 'enter' fragment for use in scripts like open_chip.
        with open(os.path.join(self.run_dir, "enter"), "w") as f:
            f.write("""
export ICC_HOME="{icc_home}"
export PATH="{icc_dir}:$PATH"
export MGLS_LICENSE_FILE="{mgls}"
export SNPSLMD_LICENSE_FILE="{snps}"
""".format(
            icc_home=dirname(dirname(icc_bin)),
            icc_dir=os.path.dirname(icc_bin),
            mgls=self.get_setting("synopsys.MGLS_LICENSE_FILE"),
            snps=self.get_setting("synopsys.SNPSLMD_LICENSE_FILE")
))

#~ if [[ "$icv" != "" ]]
#~ then
#~ cat >>"$run_dir"/enter <<EOF
#~ export ICV_HOME_DIR="$(dirname $(dirname $(dirname $icv)))"
#~ export PATH="$(dirname $icv):\$PATH"
#~ EOF
#~ fi

        # Pre-extract the Synopsys reference methodology for ICC.
        self.run_executable([
            "tar", "-xf", synopsys_rm_tarball, "-C", self.run_dir, "--strip-components=1"
        ])

        # Make sure that generated-scripts exists.
        generated_scripts_dir = os.path.join(self.run_dir, "generated-scripts")
        os.makedirs(generated_scripts_dir, exist_ok=True)

#~ if [[ "$icv" == "" ]]
#~ then
    #~ echo "ICV is not specified" >&2
    #~ # exit 1
#~ fi


#~ # Copy over any resource files if needed.
#~ resources_icc_path="$($get_config $config_db -n "" par.icc.resources_icc_path)"
#~ resources_tech_path="$($get_config $config_db -n "" par.icc.resources_tech_path)"
#~ resources_project_path="$($get_config $config_db -n "" par.icc.resources_project_path)"
#~ if [[ ! -z "$resources_icc_path" ]]; then
    #~ cp -ra "${resources_icc_path}"/* $run_dir
#~ fi
#~ if [[ ! -z "$resources_tech_path" ]]; then
    #~ cp -ra "${resources_tech_path}"/* $run_dir
#~ fi
#~ if [[ ! -z "$resources_project_path" ]]; then
    #~ cp -ra "${resources_project_path}"/* $run_dir
#~ fi

        # Gather/load libraries.
        timing_dbs = ' '.join(self.technology.read_libs([hammer_tech.filters.timing_db_filter], HammerTechnologyUtils.to_plain_item))
        milkyway_lib_dirs = ' '.join(self.technology.read_libs([hammer_tech.filters.milkyway_lib_dir_filter], HammerTechnologyUtils.to_plain_item))
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
        verilog_args = list(self.input_files)
        error = False
        for v in verilog_args:
            if not (v.endswith(".v")):
                self.logger.error("Non-verilog input {0} detected! ICC only supports Verilog for its netlists.".format(v))
                error = True
            if not os.path.isfile(v):
                self.logger.error("Input file {0} does not exist!".format(v))
                error = True
        if error:
            return False

        if len(verilog_args) < 1:
            self.logger.error("No post-synthesis Verilog netlists specified!")
            return False

        verilogs = ' '.join(verilog_args)

        # Start customizing the reference methodology.
        # Set all the input files, etc.
        with open(os.path.join(self.run_dir, "rm_setup", "common_setup.tcl"), "a") as f:
            f.write("""
set DESIGN_NAME "{top_module}";
set RTL_SOURCE_FILES "{verilogs}";
set TARGET_LIBRARY_FILES "{timing_dbs}";
set MW_REFERENCE_LIB_DIRS "{milkyway_lib_dirs}";
set MIN_LIBRARY_FILES "";
set TECH_FILE "{milkyway_techfiles}";
set MAP_FILE "{tlu_map}";
set TLUPLUS_MAX_FILE "{tlu_max_caps}";
set TLUPLUS_MIN_FILE "{tlu_min_caps}";
set ALIB_DIR "alib";
set DCRM_CONSTRAINTS_INPUT_FILE "generated-scripts/constraints.tcl";
set REPORTS_DIR "reports";
set RESULTS_DIR "results";
set CLOCK_UNCERTAINTY "0.04";
set INPUT_DELAY "0.10";
set OUTPUT_DELAY "0.10";
set ICC_NUM_CORES {num_cores};
set_host_options -max_cores {num_cores};

set MW_POWER_NET                "{MW_POWER_NET}";
set MW_POWER_PORT               "{MW_POWER_PORT}";
set MW_GROUND_NET               "{MW_GROUND_NET}";
set MW_GROUND_PORT              "{MW_GROUND_PORT}";
""".format(
                top_module=self.top_module,
                verilogs=verilogs,
                timing_dbs=timing_dbs,
                milkyway_lib_dirs=milkyway_lib_dirs,
                num_cores=self.get_setting("vlsi.core.max_threads"),
                milkyway_techfiles=milkyway_techfiles,
                tlu_map=tlu_map,
                tlu_max_caps=tlu_max_caps,
                tlu_min_caps=tlu_min_caps,
                MW_POWER_NET=self.get_setting("par.icc.MW_POWER_NET"),
                MW_POWER_PORT=self.get_setting("par.icc.MW_POWER_PORT"),
                MW_GROUND_NET=self.get_setting("par.icc.MW_GROUND_NET"),
                MW_GROUND_PORT=self.get_setting("par.icc.MW_GROUND_PORT")
            ))

        icc_setup_path = os.path.join(self.run_dir, "rm_setup", "icc_setup.tcl")
        common_setup_path = os.path.join(self.run_dir, "rm_setup", "common_setup.tcl")

        self.append_contents_to_path(str(self.get_setting("par.icc.common_setup_appendix_tcl_contents", nullvalue="")),
                                     common_setup_path)
        self.append_contents_to_path(str(self.get_setting("par.icc.icc_setup_appendix_tcl_contents", nullvalue="")),
                                     icc_setup_path)

#~ # Read the core's configuration file to figure out what all the clocks should
#~ # look like.
#~ cat >> $run_dir/generated-scripts/constraints.tcl <<"EOF"
#~ if {![info exists generated_scripts_constraints_included]} {
#~ set generated_scripts_constraints_included 1;
#~ EOF

#~ python3 >>$run_dir/generated-scripts/constraints.tcl <<EOF
#~ import json
#~ with open("${config}") as f:
    #~ config = json.load(f)

#~ import re
#~ for clock in config["clocks"]:
    #~ clock_name = clock["name"]
    #~ clock_period = clock["period"]
    #~ par_derating = clock["par derating"]
    #~ if not re.match("[0-9]+ *[np]s", clock_period):
        #~ error
    #~ if not re.match("[0-9]+ *[np]s", par_derating):
        #~ error

    #~ if re.match("[0-9]+ *ns", clock_period):
        #~ clock_period_ns = re.sub(" *ns", "", clock_period)
    #~ if re.match("[0-9]+ *ps", clock_period):
        #~ clock_period_ns = int(re.sub(" *ps", "", clock_period)) / 1000.0

    #~ if re.match("[0-9]+ *ns", par_derating):
        #~ par_derating_ns = re.sub(" *ns", "", par_derating)
    #~ if re.match("[0-9]+ *ps", par_derating):
        #~ par_derating_ns = int(re.sub(" *ps", "", par_derating)) / 1000.0

    #~ print("create_clock {0} -name {0} -period {1}".format(clock_name, clock_period_ns + par_derating_ns))
    #~ print("set_clock_uncertainty 0.01 [get_clocks {0}]".format(clock_name))
#~ EOF

#~ # The constraints file determines how the IO is constrained and what the clocks
#~ # look like.
#~ cat >> $run_dir/generated-scripts/constraints.tcl <<"EOF"
#~ # set drive strength for inputs
#~ #set_driving_cell -lib_cell INVD0BWP12T [all_inputs]
#~ # set load capacitance of outputs
#~ set_load -pin_load 0.004 [all_outputs]

#~ #set all_inputs_but_clock [remove_from_collection [all_inputs] [get_ports clock]]
#~ #set_input_delay 0.02 -clock [get_clocks clock] $all_inputs_but_clock
#~ #set_output_delay 0.03 -clock [get_clocks clock] [all_outputs]

#~ #set_isolate_ports [all_outputs] -type buffer
#~ #set_isolate_ports [remove_from_collection [all_inputs] clock] -type buffer -force
#~ EOF

#~ # We allow users to specify metal routing directions since some technologies
#~ # don't support those.
#~ python3 >>$run_dir/generated-scripts/constraints.tcl <<EOF
#~ import json
#~ with open("${technology}") as f:
    #~ config = json.load(f)

#~ # Suppress PSYN-882 ("Warning: Consecutive metal layers have the same preferred routing direction") while the layer routing is being built.
#~ print("set suppress_errors  [concat \$suppress_errors  [list PSYN-882]]")

#~ for library in config["libraries"]:
    #~ if "metal layers" in library:
        #~ for layer in library["metal layers"]:
            #~ print("set_preferred_routing_direction -layers {{ {0} }} -direction {1}".format(layer["name"], layer["preferred routing direction"]))

#~ print("set suppress_errors  [lminus \$suppress_errors  [list PSYN-882]]")
#~ EOF

#~ cat >> $run_dir/generated-scripts/constraints.tcl <<"EOF"
#~ }
#~ # generated_scripts_constraints_included
#~ EOF
        mcmm_script = self.get_setting("par.icc.mcmm_script")
        # Use DC's Verilog output instead of the milkyway stuff, which requires
        # some changes to the RM.
        # FIXME: There's a hidden dependency on the SDC file here.
        self.replace_tcl_set("ICC_INIT_DESIGN_INPUT", "VERILOG", icc_setup_path)
        self.replace_tcl_set("ICC_IN_VERILOG_NETLIST_FILE", verilogs, icc_setup_path)
        self.replace_tcl_set("ICC_MCMM_SCENARIOS_FILE", mcmm_script, icc_setup_path)
        # self.replace_tcl_set("ICC_IN_SDC_FILE", "../syn-rundir/results/{top_module}.mapped.sdc" , icc_setup_path.format(top_module = self.top_module))
        self.replace_tcl_set("ICC_NUM_CORES", self.get_setting("vlsi.core.max_threads"), icc_setup_path, quotes=False)

#~ # If there's no ICV then don't run any DRC stuff at all.
#~ if [[ "$icv" != "" ]]
#~ then
    #~ # ICC claims this is only necessary for 45nm and below, but I figure if anyone
    #~ # provides ICV metal fill rules then we might as well go ahead and use th
    #~ # rather than ICC's built-in metal filling.
    #~ if [[ "$metal_fill_ruleset" != "" ]]
    #~ then
        #~ sed 's@^set ADD_METAL_FILL.*@set ADD_METAL_FILL "ICV";@' -i $run_dir/rm_setup/icc_setup.tcl
        #~ sed "s@^set SIGNOFF_FILL_RUNSET .*@set SIGNOFF_FILL_RUNSET \"${metal_fill_ruleset[@]}\";@" -i $run_dir/rm_setup/icc_setup.tcl
    #~ fi

    #~ if [[ "$signoff_ruleset" != "" ]]
    #~ then
        #~ sed "s@^set SIGNOFF_DRC_RUNSET .*@set SIGNOFF_DRC_RUNSET \"${signoff_ruleset[@]}\";@" -i $run_dir/rm_setup/icc_setup.tcl
    #~ fi
#~ else
    #~ sed 's@^set ADD_METAL_FILL.*@set ADD_METAL_FILL "";@' -i $run_dir/rm_setup/icc_setup.tcl
#~ fi

#~ # The technology is expected to provide a list of filler cells that ICC uses.
#~ filler_metal_cells_list=$($list_macros -l $technology_macro_library -t "metal filler" | xargs echo)
#~ filler_cells_list=$($list_macros -l $technology_macro_library -t filler | xargs echo)
#~ sed 's@^set ADD_FILLER_CELL .*@set ADD_FILLER_CELL TRUE@' -i $run_dir/rm_setup/icc_setup.tcl
#~ sed "s@^set FILLER_CELL_METAL .*@set FILLER_CELL_METAL \"${filler_metal_cells_list}\";@" -i $run_dir/rm_setup/icc_setup.tcl
#~ sed "s@^set FILLER_CELL .*@set FILLER_CELL \"${filler_cells_list}\";@" -i $run_dir/rm_setup/icc_setup.tcl

#~ # I want ICC to run all the sanity checks it can
#~ sed "s@^set ICC_SANITY_CHECK.*@set ICC_SANITY_CHECK TRUE@" -i $run_dir/rm_setup/icc_setup.tcl

#~ # If I don't ask ICC for high-effort place/route then it does a bad job, so
#~ # just always ask for high effort.
#~ # FIXME: This should probably be a user tunable...
#~ sed 's@^set PLACE_OPT_EFFORT.*@set PLACE_OPT_EFFORT "high"@' -i $run_dir/rm_setup/icc_setup.tcl
#~ sed 's@^set ROUTE_OPT_EFFORT.*@set ROUTE_OPT_EFFORT "high"@' -i $run_dir/rm_setup/icc_setup.tcl
#~ sed 's@^set ICC_TNS_EFFORT_PREROUTE.*@set ICC_TNS_EFFORT_PREROUTE "HIGH"@' -i $run_dir/rm_setup/icc_setup.tcl
#~ sed 's@^set ICC_TNS_EFFORT_POSTROUTE.*@set ICC_TNS_EFFORT_POSTROUTE "HIGH"@' -i $run_dir/rm_setup/icc_setup.tcl

#~ # Some venders need the extra layer IDs.  While this is vendor-specific, I
#~ # don't see a reason to turn it off for everyone else.
#~ sed 's@^set MW_EXTENDED_LAYER_MODE.*@set MW_EXTENDED_LAYER_MODE TRUE@' -i $run_dir/rm_setup/icc_setup.tcl

#~ # For some reason, ICC isn't echoing some of my user input files.  I want it to.
#~ sed 's@source \$@source -echo $@g' -i $run_dir/rm_*/*.tcl

#~ # The only difference between this script and the actual ICC run is that this
#~ # one generates a list of macros that will be used to floorplan the design, while
#~ # the other one actually
#~ cat > $run_dir/generated-scripts/list_macros.tcl <<EOF
#~ source rm_setup/icc_setup.tcl
#~ open_mw_lib ${top}_LIB
#~ open_mw_cel -readonly init_design_icc

#~ create_floorplan -control_type aspect_ratio -core_aspect_ratio 1 -core_utilization 0.7 -left_io2core 3 -bottom_io2core 3 -right_io2core 3 -top_io2core 3 -start_first_row

#~ set top_left_x [lindex [get_placement_area] 0]
#~ set top_left_y [lindex [get_placement_area] 1]
#~ set bottom_right_x [lindex [get_placement_area] 2]
#~ set bottom_right_y [lindex [get_placement_area] 3]
#~ echo "${top} module=${top} top_left=(\$top_left_x, \$top_left_y) bottom_right=(\$bottom_right_x, \$bottom_right_y)" >> results/${top}.macros.out

#~ set fixed_cells [get_fp_cells -filter "is_fixed == true"]
#~ foreach_in_collection cell \$fixed_cells {
    #~ set full_name [get_attribute \$cell full_name]
    #~ set ref_name [get_attribute \$cell ref_name]
    #~ set height [get_attribute \$cell height]
    #~ set width [get_attribute \$cell width]
    #~ echo "\$full_name parent=${top} module=\$ref_name width=\$width height=\$height" >> results/${top}.macros.out
#~ }
#~ exit
#~ EOF

#~ ## FIXME: This throws errors because it's accessing some views on disk.
#~ ## I want ICC to try and fix DRCs automatically when possible.  Most of the
#~ ## commands are commented out for some reason, this enables them.
#~ #drc_runset="$(echo "${signoff_ruleset[@]}" | xargs basename -s .rs)"
#~ #sed 's@^set ICC_ECO_SIGNOFF_DRC_MODE .*@set ICC_ECO_SIGNOFF_DRC_MODE "AUTO_ECO"@' -i $run_dir/rm_setup/icc_setup.tcl
#~ #sed 's@.*#  s\(.*\)@s\1@' -i $run_dir/rm_icc_zrt_scripts/signoff_drc_icc.tcl
#~ #sed "s@^\\(signoff_autofix_drc .*\\)@exec $ICV_HOME_DIR/contrib/generate_layer_rule_map.pl -dplog signoff_drc_run/run_details/$drc_runset.dp.log -tech_file $(readlink -f ${tf[@]}) -o signoff_autofix_drc.config\n\\1@" -i $run_dir/rm_icc_zrt_scripts/signoff_drc_icc.tcl
#~ #sed 's@\$config_file@signoff_autofix_drc.config@' -i $run_dir/rm_icc_zrt_scripts/signoff_drc_icc.tcl

#~ # FIXME: I actually can't insert double vias on SAED32 becaues of DRC errors.
#~ # It smells like the standard cells just aren't setup for it, but this needs to
#~ # be fixed somehow as it'll be necessary for a real chip to come back working.
#~ sed 's@set ICC_DBL_VIA .*@set ICC_DBL_VIA FALSE@' -i $run_dir/rm_setup/icc_setup.tcl
#~ sed 's@set ICC_DBL_VIA_FLOW_EFFORT .*@set ICC_DBL_VIA_FLOW_EFFORT "NONE"@' -i $run_dir/rm_setup/icc_setup.tcl
        
       
        # Floorplan Generation
        floorplan_mode = self.get_setting("par.icc.floorplan_mode")
        floorplan_script = self.get_setting("par.icc.floorplan_script")

        if floorplan_mode == "annotations":
            raise NotImplementedError("Not implemented")
        elif floorplan_mode == "manual":
            if floorplan_script == "null" or floorplan_script == "":
                self.logger.error("floorplan_mode is manual but no floorplan_script specified")
                return False
            copyfile(floorplan_script, os.path.join(self.run_dir, "generated-scripts", "floorplan_inner.tcl"))
            self.replace_tcl_set("ICC_FLOORPLAN_INPUT", "USER_FILE", icc_setup_path)
            self.replace_tcl_set("ICC_IN_FLOORPLAN_USER_FILE", "generated-scripts/floorplan_inner.tcl", icc_setup_path)
        elif floorplan_mode == "generate":
            floorplan_script = os.path.join(self.run_dir, "generated-scripts", "floorplan_inner.tcl")
            with open(floorplan_script, "w") as f:
                f.write("\n".join(self.generate_floorplan_tcl()))
            self.replace_tcl_set("ICC_FLOORPLAN_INPUT", "USER_FILE", icc_setup_path)
            self.replace_tcl_set("ICC_IN_FLOORPLAN_USER_FILE", "generated-scripts/floorplan_inner.tcl", icc_setup_path)
        elif floorplan_mode == "default":
            self.replace_tcl_set("ICC_FLOORPLAN_INPUT", "CREATE", icc_setup_path)
        elif floorplan_mode == "auto":
            ## TODO : need to check if this mode is feasible depending on whether ICC has a command similar to "plan_design" in Innovus
            pass        
        else:
            self.logger.error("Invalid floorplan_mode %s" % (floorplan_mode))
            return False

        # Prepend constraints to the floorplan.
        with open(os.path.join(self.run_dir, "generated-scripts", "floorplan.tcl"), "w") as f:
            f.write("""
source -echo -verbose generated-scripts/constraints.tcl
source -echo -verbose generated-scripts/floorplan_inner.tcl
""")
        
        icc_init_design = os.path.join(self.run_dir, "rm_icc_scripts", "init_design_icc.tcl")
        
        # Bumps Placement
        bumps_mode = self.get_setting("vlsi.inputs.bumps_mode")
        if bumps_mode == "manual":
            bumps = self.get_bumps()
            if bumps is not None:  
                #Call the Bump-Placement API method
                with open(icc_init_design, "a") as f:
                        f.write("\n".join(self.place_bumps()))
            else:
                self.logger.error("Bumps Placement Mode is manual but no bumps are specified!!")
                return False
        elif bumps_mode == "empty":
            self.logger.warning("Bumps Placement Mode is empty. No bumps are added.")
        else:
            self.logger.error("Invalid Bumps Placement Mode %s" % (bumps_mode))
            return False
 
        #TapCell Placement
        tap_cells = self.technology.get_special_cell_by_type(CellType.TapCell)
        if len(tap_cells) == 0:
            self.logger.warning("Tap cells are improperly defined in the tech plugin and will not be added."
                                "This step should be overridden with a user hook.")
        else:
            #Call the Tap-Cell Placement API method  
            with open(icc_init_design, "a") as f:
                    f.write("\n".join(self.place_tap_cells()))
        
        # Pins Placement
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
                    with open(icc_init_design, "a") as f:
                            f.write("\n".join(self.place_pins()))
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
      
        # Power Straps Generation
        power_straps_mode = self.get_setting("par.power_straps_mode")
        if power_straps_mode == "empty":
            self.logger.warning("Power Straps Mode is empty!! No power straps are generated")
        elif power_straps_mode == "manual":
            # Reads the manual power straps script contents and appends it to icc_init_design
            power_straps_script = self.get_setting("par.power_straps_script_contents")
            power_script_iterable = [] # List[str] - list containing the individual lines of the input power straps script
            with open(power_straps_script, "r") as f:
                    for line in f:
                        stripped_line = line.strip()
                        power_script_iterable.append(stripped_line)
            # Appending the input power script to icc_init_design file
            with open(icc_init_design, "a") as f:
                    f.write("\n".join(power_script_iterable))     
        elif power_straps_mode == "generate":
            # Generate the Power Straps by passing the input from Hammer IR to the Power Straps Generation API
            # Call the Power Straps API methods
            with open(icc_init_design, "a") as f:
                    f.write("\n".join(self.create_power_straps_tcl()))
        else:
            self.logger.error("Invalid Power Straps Mode %s" % (power_straps_mode))
            return False
          
#~ # Opens the floorplan straight away, which is easier than doing it manually
#~ cat > $run_dir/generated-scripts/open_floorplan.tcl <<EOF
#~ source rm_setup/icc_setup.tcl
#~ open_mw_lib -r ${top}_LIB
#~ open_mw_cel -r floorplan
#~ EOF

#~ cat > $run_dir/generated-scripts/open_floorplan <<EOF
#~ cd $run_dir
#~ source enter
#~ $ICC_HOME/bin/icc_shell -gui -f generated-scripts/open_floorplan.tcl
#~ EOF
#~ chmod +x $run_dir/generated-scripts/open_floorplan

#~ cat > $run_dir/generated-scripts/open_power.tcl <<EOF
#~ source rm_setup/icc_setup.tcl
#~ open_mw_lib -r ${top}_LIB
#~ open_mw_cel -r power
#~ EOF

#~ cat > $run_dir/generated-scripts/open_power <<EOF
#~ cd $run_dir
#~ source enter
#~ $ICC_HOME/bin/icc_shell -gui -f generated-scripts/open_power.tcl
#~ EOF
#~ chmod +x $run_dir/generated-scripts/open_power

        with open(os.path.join(generated_scripts_dir, "open_chip.tcl"), "w") as f:
            f.write("""
cd {run_dir}
source rm_setup/icc_setup.tcl
open_mw_lib -r {top}_LIB
open_mw_cel -r chip_finish_icc
""".format(run_dir=self.run_dir, top=self.top_module))

        with open(os.path.join(generated_scripts_dir, "open_chip"), "w") as f:
            f.write("""
cd {run_dir}
source enter
$ICC_HOME/bin/icc_shell -gui -f generated-scripts/open_chip.tcl
""".format(run_dir=self.run_dir))
        self.run_executable([
            "chmod", "+x", os.path.join(generated_scripts_dir, "open_chip")
        ])

#~ # Write SDF
#~ cat > $run_dir/generated-scripts/write_sdf.tcl <<EOF
#~ source rm_setup/icc_setup.tcl
#~ open_mw_cel \$ICC_OUTPUTS_CEL -lib \$MW_DESIGN_LIBRARY
#~ current_design ${top}
#~ write_sdf \$RESULTS_DIR/\$DESIGN_NAME.output.sdf
#~ exit
#~ EOF

#~ cat > $run_dir/generated-scripts/write_sdf <<EOF
#~ cd $run_dir
#~ $ICC_HOME/bin/icc_shell -f generated-scripts/write_sdf.tcl
#~ EOF
#~ chmod +x $run_dir/generated-scripts/write_sdf

        # Build args.
        args = [
            os.path.join(self.tool_dir, "tools", "run-par"),
            "--run_dir", self.run_dir
            #~ "--dc", dc_bin,
            #~ "--MGLS_LICENSE_FILE", self.get_setting("synopsys.MGLS_LICENSE_FILE"),
            #~ "--SNPSLMD_LICENSE_FILE", self.get_setting("synopsys.SNPSLMD_LICENSE_FILE"),
            #~ "--preferred_routing_directions_fragment", preferred_routing_directions_fragment,
            #~ "--find_regs_tcl", os.path.join(self.tool_dir, "tools", "find-regs.tcl"),
            #~ "--top", self.top_module
        ]

        # Temporarily disable colours/tag to make DC run output more readable.
        # TODO: think of a more elegant way to do this?
        HammerVLSILogging.enable_colour = False
        HammerVLSILogging.enable_tag = False
        self.run_executable(args, cwd=self.run_dir) # TODO: check for errors and deal with them
        HammerVLSILogging.enable_colour = True
        HammerVLSILogging.enable_tag = True

        return True


tool = ICC
