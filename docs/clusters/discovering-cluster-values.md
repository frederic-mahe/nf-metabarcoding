# Discovering the values for a cluster profile

`[S87]`. This is the recipe for filling in (or checking) a
`conf/clusters/<name>.config` from a real slurm cluster. Run the commands
below **on a login node** of the target cluster, then map each result to
the knob it sets (table at the end).

Why bother instead of copying an nf-core config? Because they go stale.
Every nf-core-derived value in this repo's French profiles was wrong when
checked against the live scheduler — `abims`, `genotoul` and `ifb_core`
were each off by 2–10× on cores and/or memory, and `genotoul`'s
container engine was reported as absent when it is actually a module. A
five-minute check on the hardware is the difference between a profile
that schedules and one that silently mis-sizes every job.

> Confirm on the hardware; do not trust a config (ours included) you have
> not checked. Partition policies, node sizes and module names all drift.


## 1. Partitions, node sizes, wall-time limits

```bash
# partition -> max/default time, cores, mem(MB), GPUs, node count
sinfo -o "%P %.14l %.14L %.6c %.10m %.10G %.6D"
#   %P partition  %l MaxTime  %L DefaultTime  %c cores/node
#   %m mem/node(MB)  %G GRES(gpus)  %D node count

# authoritative per-partition limits, incl. AllowAccounts / MaxTime
scontrol show partition

# the single largest node (sort by mem, then cores) -> the ceiling
sinfo -N -o "%.20N %.6c %.10m %.16P" | sort -k3 -n | tail
```

What to read off:

- **The biggest node** (most memory, most cores) sets `max_memory` /
  `max_cpus` / `resourceLimits`. Set `max_memory` *just under* the node's
  MB (e.g. a 2064000 MB node → `2000.GB`) so a request can actually
  schedule — RAM in MB ÷ 1024 = GiB, and Nextflow's `GB`/`TB` are binary
  (GiB/TiB), matching how slurm reports memory.
- **A node can belong to several partitions.** Check the `Nodes=` list:
  if the big-memory node is already in the default partition (as on
  `genotoul`, where `bigmem01` is in `workq`), you do **not** need a
  separate memory route. If the big-memory partition holds genuinely more
  RAM than the default (as on `abims`/`saga`), route to it by
  `task.memory` in the `queue` closure.
- **`MaxTime` per partition** drives the `queue` time cutoff (e.g.
  `workq` MaxTime 4 d → `task.time <= 96.h ? 'workq' : 'unlimitq'`).
- **`DefaultTime`** is informational: this pipeline sets an explicit
  `time` on every process tier (`conf/slurm.config`), so a low site
  default (e.g. genotoul's 2 h) never bites.
- **`AllowAccounts`**: `ALL` means `--slurm_account` is optional; a
  specific list (or a per-partition account, as on `meso`) means it is
  required — see step 2.
- Ignore GPU / docker / interactive / dev partitions (this pipeline has
  no GPU steps).


## 2. Accounts you may charge to

```bash
sacctmgr show assoc user=$USER format=Account%25,Partition%20,QOS%20,MaxWall,GrpTRES%30
sacctmgr show user $USER withassoc format=User,Account%25,Partition%20
sshare -U          # fairshare accounts you belong to
```

- The **account string(s)** you are allowed to use → pass with
  `--slurm_account` (most sites), or hardcode in the profile only when
  the site *derives* the account from the job (as `meso` does: bigmem vs
  cpu, and a wall-time tier).
- The **per-account `MaxWall`** can encode the policy: on `meso`,
  `dedicated-cpu@cirad-default` = 1 h, `-normal` = 2 d, `-long` = 120 d —
  which is exactly the account-routing closure in
  `conf/clusters/meso.config`.


## 3. Scratch / temporary storage

These variables only exist *inside* an allocation, so probe with a tiny
job:

```bash
srun --time=2:00 --mem=1G --pty bash -c \
  'echo "TMPDIR=$TMPDIR LOCALSCRATCH=$LOCALSCRATCH SCRATCH=$SCRATCH USERWORK=$USERWORK"'
```

Whichever is populated and points at fast per-job storage can become
`process.scratch` in the profile (e.g. saga's `$SCRATCH`). If only
`$TMPDIR` is set (as on `ifb_core`/`abims`), leave `scratch` unset — tasks
run in the work dir on the shared filesystem. Note `$LOCALSCRATCH` often
needs an explicit `--gres=localscratch:<size>` request, so it is not a
free default.


## 4. Container engine + bind mounts

```bash
# IMPORTANT: engines may be on PATH *or* hidden behind a namespaced
# module. Check both, and search the whole module tree:
which -a singularity apptainer 2>/dev/null
module -t avail 2>&1 | grep -iE 'apptainer|singularity|containers/'

# the filesystems worth bind-mounting (the roots your data/work live on)
findmnt -lo TARGET,FSTYPE | grep -iE 'scratch|work|shared|cluster|bank|save|projet|home'
```

- **On PATH** (`/usr/bin/singularity`, as on abims/ifb_core) → just set
  `singularity`/`apptainer` `autoMounts` + `runOptions` binds in the
  profile.
- **Module only** (as on `genotoul`: `containers/Apptainer/1.4.1`,
  `containers/singularity/CE-4.3.1`) → load it via an engine-aware
  `process.beforeScript` so the binary is on PATH when Nextflow launches
  the container. See `conf/clusters/genotoul.config` for the pattern.
  > This is the step the first-pass commands missed: `module avail
  > singularity` finds nothing because the module is under the
  > `containers/` namespace, and `command -v` fails because it is not
  > loaded. Always grep the full `module avail` tree.
- **Neither** → use `-profile <name>,conda` or `-profile <name>,modules`.
- **Bind mounts**: bind the shared roots from `findmnt` (e.g. `-B /shared`,
  `-B /work -B /save`). `autoMounts` covers the work dir; explicit binds
  are for inputs/references on other roots.
- **Image cache**: point `cacheDir` at a **writable** area you own
  (a project/scratch dir), never an admin-managed software tree.


## 5. Tool modules (for `-profile <name>,modules`)

```bash
module -t avail 2>&1 | grep -iE 'vsearch|swarm|cutadapt|mumu'
```

Feeds `--module_vsearch` / `--module_swarm` / `--module_cutadapt` /
`--module_mumu`. Watch for gaps: `genotoul` has no `mumu` module, so a
container (where `mumu` comes from `environment.yml` via Wave) or conda is
simpler there.


## Value → knob cheat-sheet

| Discovered from | Sets in `conf/clusters/<name>.config` |
|-----------------|----------------------------------------|
| largest node mem / cores (`sinfo`) | `params.max_memory`, `params.max_cpus`, `resourceLimits` |
| partition `MaxTime` + which holds big-mem nodes | `process.queue` closure (time and/or `task.memory`) |
| `AllowAccounts` / `sacctmgr` accounts + `MaxWall` | `--slurm_account`, or the `clusterOptions` account closure (meso) |
| scratch env var (`srun` probe) | `process.scratch` |
| engine on PATH vs module (`which`/`module avail`) | `singularity`/`apptainer` blocks, or `process.beforeScript` module-load |
| shared FS roots (`findmnt`) | `singularity.runOptions` / `apptainer.runOptions` binds |
| tool module names | `--module_vsearch` / `_swarm` / `_cutadapt` / `_mumu` |

Once filled in, register the profile in `nextflow.config` (two
`includeConfig` lines — see `conf/clusters/_template.config`) and add a
row to `tests/check-cluster-profiles.sh` so `nextflow config -profile
<name>` is checked in CI.
