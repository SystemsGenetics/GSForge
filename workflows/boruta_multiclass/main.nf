#!/usr/bin/env nextflow


// Create a helper function that will produce argument combinations.
import groovy.json.JsonOutput

def create_kwargs_combinations( Map kwarg_map ) {
    kwarg_map.values().combinations { args ->
        [kwarg_map.keySet().asList(), args].transpose().collectEntries { [(it[0]): it[1]] }
    }
}

println """\

================================
Boruta Multiclass Gene Selection
================================

Workflow Information:
---------------------
  Project Directory:  ${workflow.projectDir}
  Launch Directory:   ${workflow.launchDir}
  Work Directory:     ${workflow.workDir}
  Config Files:       ${workflow.configFiles}
  Container Engine:   ${workflow.containerEngine}
  Profile(s):         ${workflow.profile}

"""

// Generate boruta argument combinations.
boruta_opts_ch = Channel.from( create_kwargs_combinations( params.boruta_opts ) )
model_opts_ch = Channel.from( create_kwargs_combinations( params.ranking_model_opts ) )


// The feature extraction process.
process boruta {

    tag { "${y_label}-${repeat}" }
    publishDir params.out_dir , mode: "move"

    input:
        each repeat from 1..params.repeats
        file( gem_netcdf ) from Channel.fromPath( params.gem_netcdf )
        each boruta_opts from boruta_opts_ch
        each model_opts from model_opts_ch
        each y_label from params.y_label
        val x_label from params.x_label
        val model from params.ranking_model

    output:
        file ('*.nc') optional true

    script:
    """
    boruta_multiclass.py --gem_netcdf ${gem_netcdf} \
                  --x_label ${x_label} \
                  --y_label ${y_label} \
                  --boruta_opts '${JsonOutput.toJson(boruta_opts)}' \
                  --ranking_model ${model} \
                  --ranking_model_opts '${JsonOutput.toJson(model_opts)}' \
                  --append_wd
    """
}
