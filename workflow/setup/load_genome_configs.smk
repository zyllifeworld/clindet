##
## Load per-genome configuration from workflow/config/genomes/.
##
## Each <genome_version>.yaml provides:
##     resources: reference genome + databases for that genome_version
##     tools:     genome-specific tool parameters/databases
## _shared.yaml holds cross-genome entries (currently `varanno`).
##
## The loader rebuilds exactly the structure the rules already use:
##     config['resources'][genome_version]
##     config['softwares_params'][genome_version]
##
## To switch over, remove these lines from wrapper/{wes,wgs,rna}.smk:
##     configfile: "workflow/config/conf/genomes.yaml"
##     configfile: "workflow/config/conf/softwares_params.yaml"
## and remove the softwares_params.yaml configfile from
## workflow/setup/rules/env_setup.smk, then include this file from the main
## Snakefile (or from a wrapper as "../workflow/setup/load_genome_configs.smk").
##
import os
import yaml
from pathlib import Path

GENOME_CONFIG_DIR = Path(workflow.basedir) / "workflow" / "config" / "genomes"
RESOURCE_ROOT = os.environ.get("CLINDET_RESOURCE_ROOT", "")


def _expand_resource_root(value):
    if isinstance(value, str):
        return value.replace("${RESOURCE_ROOT}", RESOURCE_ROOT)
    if isinstance(value, dict):
        return {k: _expand_resource_root(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand_resource_root(v) for v in value]
    return value


if not GENOME_CONFIG_DIR.is_dir():
    raise WorkflowError(
        f"Genome config directory not found: {GENOME_CONFIG_DIR}. "
        "Run workflow/setup/scripts/split_genome_configs.py first."
    )

_shared_file = GENOME_CONFIG_DIR / "_shared.yaml"
if _shared_file.exists():
    _shared = _expand_resource_root(yaml.safe_load(_shared_file.read_text()) or {})
    for _key, _value in _shared.get("resources", {}).items():
        config.setdefault("resources", {}).setdefault(_key, _value)

for _file in sorted(GENOME_CONFIG_DIR.glob("*.yaml")):
    if _file.name.startswith("_"):
        continue
    _doc = _expand_resource_root(yaml.safe_load(_file.read_text()) or {})
    _gv = _doc["genome_version"]
    config.setdefault("resources", {})[_gv] = _doc.get("resources", {})
    config.setdefault("softwares_params", {})[_gv] = _doc.get("tools", {})
