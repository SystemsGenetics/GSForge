#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Handles a single run of Boruta_Py feature selection.

This script is designed to be run with Nextflow, since each run requires
some time and a large amount of RAM. It models how to use Nextflow to
explore parameters for boruta or a ranking model.


Testing the script:

```
boruta_multiclass.py -g oryza_sativa.nc -y Treatment
```

Testing with nextflow:

```
nextflow run boruta_multiclass.nf
```

# TODO:
+ Get lgb models working.

"""

import os
import click
import xarray as xr
import json
import hashlib

from GSForge.models._AnnotatedGEM import AnnotatedGEM  # To avoid the __init__ run.
from GSForge.operations.prospectors import boruta_prospector

from sklearn.ensemble import RandomForestClassifier
# import lightgbm as lgb


MODELS = {
    "RandomForestClassifier": RandomForestClassifier,
    # "LGBMClassifier": lgb.LGBMClassifier,
}

MODEL_DEFAULTS = {
    "RandomForestClassifier": '{"class_weight": "balanced", "max_depth": 3, "n_jobs: -1}',
    # "LGBMClassifier": 'null',
}


@click.command()
@click.option('-g', '--gem_netcdf', 'gem_netcdf')
@click.option('-x', '--x_label', 'x_label', default='counts')
@click.option('-y', '--y_label', 'y_label')
@click.option('-bo', '--boruta_opts', 'boruta_opts', default='{"perc": 100, "max_iter": 2}')
@click.option('-rm', '--ranking_model', 'ranking_model', default='RandomForestClassifier')
@click.option('-rmo', '--ranking_model_opts', 'ranking_model_opts', default=None)
@click.option('--append_wd/--no-append_wd', default=False)
def main(gem_netcdf,
         x_label,
         y_label,
         boruta_opts,
         ranking_model_opts,
         ranking_model,
         append_wd,
         ):
    """Run standard boruta_py feature selection with one of the following models:
        * RandomForestClassifier
        * LGBMClassifier
    """
    click.echo("Boruta runner started.")

    # Construct a hash of all given arguments.
    def h6(w):
        h = hashlib.md5(w.encode())
        return h.hexdigest()[:6]

    arg_str = "_".join([
        gem_netcdf,
        x_label,
        y_label,
        boruta_opts,
        ranking_model_opts,
        ranking_model,
    ])
    arg_hash = h6(arg_str)

    click.echo("Parsing options...")
    if ranking_model_opts is None:
        ranking_model_opts = MODEL_DEFAULTS[ranking_model]

    model_kwargs = json.loads(ranking_model_opts)
    boruta_kwargs = json.loads(boruta_opts)

    click.echo("Constructing Model...")
    selection_model = MODELS[ranking_model](**model_kwargs)

    click.echo("Loading Data...")
    # Load the data in a safe way so that multiple processes can read the file.
    with xr.open_dataset(gem_netcdf) as ds:
        ds.load()
        agem = AnnotatedGEM(data=ds)

        click.echo("Running boruta feature selection...")
        boruta_results = boruta_prospector(
            agem,
            estimator=selection_model,
            annotation_variables=y_label,
            **boruta_kwargs,)
        # Add the argument hash to the attributes.
        boruta_results = boruta_results.assign_attrs({"arg_hash": arg_hash})

    click.echo("Examining Results...")
    found_gene_count = boruta_results.support.sum().values
    if found_gene_count > 0:
        click.echo(f"{found_gene_count} Genes found.", color="green")
    else:
        weak_gene_count = boruta_results.support_weak.sum().values
        if weak_gene_count > 0:
            click.echo(f"{weak_gene_count} genes with weak support found.", color="yellow")
        else:
            click.echo(f"No genes selected for target {y_label}.", color="red")

    if append_wd:
        output_dir = os.path.basename(os.getcwd())
        output_filename = f"{arg_hash}_{output_dir[-6:]}.nc"
    else:
        output_filename = f"{arg_hash}.nc"

    boruta_results.to_netcdf(output_filename)
    click.echo(f"Results saved to: {output_filename}")

    with mlflow.start_run():
        mlflow.log_param("x", 1)
        mlflow.log_metric("y", 2)


    exit(0)


if __name__ == '__main__':
    main()
